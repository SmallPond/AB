// Copyright (c) 2025 -  dingmos
// Use of this source code is governed by a GPLv3 license that can be found in the LICENSE file.

/*
    SCARA Kinematics - Absolute Angle Mode
    
    This implements SCARA kinematics with absolute angle mode where:
    - theta (first arm angle): measured from the -X axis
    - psi (second arm angle): measured from world coordinate system
    
    This is different from relative angle SCARA where the second arm
    angle is measured relative to the first arm.
*/

#include "Scara.h"

#include "Machine/MachineConfig.h"
#include "Limit.h"
#include "Machine/Homing.h"
#include "Protocol.h"
#include "GCode.h"    // gc_sync_position
#include "Settings.h" // coords, CoordIndex for G28 position
#include "Channel.h"  // Channel for log_msg_to

#include <cmath>

namespace Kinematics {

    /*
      Configuration group handler
      Defines the YAML configuration parameters
    */
    void Scara::group(Configuration::HandlerBase& handler) {
        handler.item("linkage1_mm", _linkage1_mm, 10.0f, 500.0f);
        handler.item("linkage2_mm", _linkage2_mm, 10.0f, 500.0f);
        handler.item("scara_offset_x", _scara_offset_x, -500.0f, 500.0f);
        handler.item("scara_offset_y", _scara_offset_y, -500.0f, 500.0f);
        handler.item("kinematic_segment_len_mm", _kinematic_segment_len_mm, 0.1f, 20.0f);
        handler.item("cal_init_theta", _cal_init_theta, -180.0f, 180.0f);
        handler.item("cal_init_psi", _cal_init_psi, -180.0f, 360.0f);
    }

    /*
      Initialize the kinematic system
    */
    void Scara::init() {
        log_info("Kinematic system: " << name());
        log_info("  Linkage 1: " << _linkage1_mm << " mm");
        log_info("  Linkage 2: " << _linkage2_mm << " mm");
        log_info("  Offset X: " << _scara_offset_x << " mm");
        log_info("  Offset Y: " << _scara_offset_y << " mm");

        init_position();
    }

    /*
      Initialize the machine position
    */
    void Scara::init_position() {
        auto  n_axis = Axes::_numberAxis;
        float min_mpos[MAX_N_AXIS];
        float max_mpos[MAX_N_AXIS];

        for (axis_t axis = X_AXIS; axis < n_axis; axis++) {
            set_steps(axis, 0);
            min_mpos[axis] = limitsMinPosition(axis);
            max_mpos[axis] = limitsMaxPosition(axis);
        }
        
        transform_cartesian_to_motors(_min_motor_pos, min_mpos);
        transform_cartesian_to_motors(_max_motor_pos, max_mpos);
        
        // Initialize last motor position
        for (axis_t axis = X_AXIS; axis < n_axis; axis++) {
            _last_motor_pos[axis] = 0.0f;
        }
    }

    /*
      Forward Kinematics: Convert motor angles to cartesian coordinates
      
      This is the "正解" algorithm from the reference:
      - theta: angle of first arm from -X axis (in degrees)
      - psi: angle of second arm from world coordinate system (in degrees)
      
      The formulas:
      x_sin = sin(theta) * L1
      x_cos = cos(theta) * L1
      y_sin = sin(psi) * L2
      y_cos = cos(psi) * L2
      
      world_x = -x_cos - y_cos
      world_y = x_sin + y_sin
      
      user_x = world_x - offset_x
      user_y = world_y - offset_y
    */
    void Scara::forward_kinematics(float* cartesian, float theta_deg, float psi_deg) {
        float theta_rad = degrees_to_radians(theta_deg);
        float psi_rad = degrees_to_radians(psi_deg);
        
        float x_sin = sinf(theta_rad) * _linkage1_mm;
        float x_cos = cosf(theta_rad) * _linkage1_mm;
        float y_sin = sinf(psi_rad) * _linkage2_mm;
        float y_cos = cosf(psi_rad) * _linkage2_mm;
        
        // World coordinates (base coordinate system)
        float world_x = -x_cos - y_cos;
        float world_y = x_sin + y_sin;
        
        // Convert to user coordinates
        cartesian[X_AXIS] = world_x - _scara_offset_x;
        cartesian[Y_AXIS] = world_y - _scara_offset_y;
    }

    /*
      Inverse Kinematics: Convert cartesian coordinates to motor angles
      
      This is the "反解" algorithm from the reference.
      
      The formulas:
      1. Convert user coords to world coords:
         world_x = -user_x - offset_x
         world_y = user_y + offset_y
      
      2. Calculate intermediate values:
         C2 = (world_x^2 + world_y^2 - L1^2 - L2^2) / (2 * L1 * L2)
         S2 = sqrt(1 - C2^2)
         K1 = L1 + L2 * C2
         K2 = L2 * S2
      
      3. Calculate angles:
         theta = atan2(K1, K2) - atan2(world_x, world_y)
         psi = atan2(S2, C2) + theta
      
      Returns true if position is reachable, false otherwise.
    */
    bool Scara::inverse_kinematics(float cartesian_x, float cartesian_y, float& theta_deg, float& psi_deg) {
        // Convert from user coordinates to world coordinates
        // Note: X is inverted in world coordinate system
        float world_x = -cartesian_x - _scara_offset_x;
        float world_y = cartesian_y + _scara_offset_y;
        
        float L1_sq = _linkage1_mm * _linkage1_mm;
        float L2_sq = _linkage2_mm * _linkage2_mm;
        
        // Distance squared from origin to target point
        float dist_sq = world_x * world_x + world_y * world_y;
        
        // Check if point is reachable
        // Must be within (L1 - L2)^2 to (L1 + L2)^2
        float min_dist = _linkage1_mm - _linkage2_mm;
        float max_dist = _linkage1_mm + _linkage2_mm;
        
        if (dist_sq < min_dist * min_dist || dist_sq > max_dist * max_dist) {
            return false;  // Point is unreachable
        }
        
        // Calculate C2 and S2 using law of cosines
        float C2 = (dist_sq - L1_sq - L2_sq) / (2.0f * _linkage1_mm * _linkage2_mm);
        
        // Clamp C2 to valid range to avoid numerical errors
        if (C2 > 1.0f) C2 = 1.0f;
        if (C2 < -1.0f) C2 = -1.0f;
        
        float S2 = sqrtf(1.0f - C2 * C2);
        
        // Calculate K1 and K2
        float K1 = _linkage1_mm + _linkage2_mm * C2;
        float K2 = _linkage2_mm * S2;
        
        // Calculate angles in radians
        float theta_rad = atan2f(K1, K2) - atan2f(world_x, world_y);
        float psi_rad = atan2f(S2, C2) + theta_rad;
        
        // Convert to degrees
        theta_deg = radians_to_degrees(theta_rad);
        psi_deg = radians_to_degrees(psi_rad);
        
        return true;
    }

    /*
      Transform cartesian coordinates to motor coordinates
      This is used for soft limit checks and position display
    */
    bool Scara::transform_cartesian_to_motors(float* motors, float* cartesian) {
        float theta_deg, psi_deg;
        
        // X and Y are transformed by inverse kinematics
        if (!inverse_kinematics(cartesian[X_AXIS], cartesian[Y_AXIS], theta_deg, psi_deg)) {
            return false;  // Position is unreachable
        }
        
        motors[X_AXIS] = theta_deg;
        motors[Y_AXIS] = psi_deg;
        
        // Other axes pass through unchanged
        auto n_axis = Axes::_numberAxis;
        for (axis_t axis = Z_AXIS; axis < n_axis; axis++) {
            motors[axis] = cartesian[axis];
        }
        
        return true;
    }

    /*
      Convert motor coordinates to cartesian coordinates
      This is used for position reporting (DRO)
    */
    void Scara::motors_to_cartesian(float* cartesian, float* motors, axis_t n_axis) {
        // X and Y are transformed by forward kinematics
        // motors[X_AXIS] is theta (degrees)
        // motors[Y_AXIS] is psi (degrees)
        forward_kinematics(cartesian, motors[X_AXIS], motors[Y_AXIS]);
        
        // Z axis is unchanged, but note that Z in cartesian comes from Z motor
        // We need to handle this separately since forward_kinematics only sets X and Y
        // Actually, forward_kinematics only touches X and Y, so we need to copy Z
        // But the base forward_kinematics already handles this correctly
        
        // Other axes pass through unchanged
        for (axis_t axis = Z_AXIS; axis < n_axis; axis++) {
            cartesian[axis] = motors[axis];
        }
    }

    /*
      cartesian_to_motors: Main motion transformation function
      
      This function is called for every linear move (G0, G1).
      It segments the cartesian move into small segments and transforms
      each segment to motor coordinates to maintain straight lines.
    */
    bool Scara::cartesian_to_motors(float* target, plan_line_data_t* pl_data, float* position) {
        auto n_axis = Axes::_numberAxis;
        
        // Save the original feed rate
        float feed_rate = pl_data->feed_rate;
        
        // Check if the destination is reachable
        float motors[MAX_N_AXIS];
        if (!transform_cartesian_to_motors(motors, target)) {
            log_warn("SCARA kinematics error. Target unreachable ("
                     << target[X_AXIS] << ", " << target[Y_AXIS] << ")");
            return false;
        }
        
        // Calculate the cartesian move distance
        float d[MAX_N_AXIS];
        copyAxes(d, target, n_axis);
        subtractAxes(d, position, n_axis);
        
        // Determine the number of segments needed for X and Y motion
        // Only X and Y need segmentation due to non-linearity
        float xyz_dist = vector_length(d, 3);
        
        uint32_t segment_count = (uint32_t)ceilf(xyz_dist / _kinematic_segment_len_mm);
        if (segment_count < 1) {
            segment_count = 1;
        }
        
        // Calculate the all-axis segment distance for feedrate scaling
        float segment_dist = vector_length(d, n_axis) / segment_count;
        
        // Calculate the cartesian segment delta
        float delta_d[MAX_N_AXIS];
        copyArray(delta_d, d, n_axis);
        multiplyArray(delta_d, 1.0f / segment_count, n_axis);
        
        // Initialize segment starting position
        float seg_position[MAX_N_AXIS];
        copyAxes(seg_position, position, n_axis);
        
        // Process each segment
        for (uint32_t segment = 1; segment <= segment_count; segment++) {
            if (sys.abort()) {
                return true;
            }
            
            // Calculate the segment end position in cartesian space
            addAxes(seg_position, delta_d, n_axis);
            
            // Transform to motor coordinates
            float seg_motors[MAX_N_AXIS];
            if (!transform_cartesian_to_motors(seg_motors, seg_position)) {
                log_error("SCARA kinematic error at position ("
                          << seg_position[X_AXIS] << ", " << seg_position[Y_AXIS] << ")");
                return false;
            }
            
            // Adjust feedrate based on motor distance vs cartesian distance
            if (!pl_data->motion.rapidMotion) {
                float motor_dist = vector_distance(seg_motors, _last_motor_pos, n_axis);
                pl_data->feed_rate = feed_rate * motor_dist / segment_dist;
            }
            
            // Send the motor move command
            if (!mc_move_motors(seg_motors, pl_data)) {
                return false;
            }
            
            // Save motor position for next segment
            copyAxes(_last_motor_pos, seg_motors, n_axis);
        }
        
        return true;
    }

    /*
      Check if a line is valid (within work envelope)
    */
    bool Scara::invalid_line(float* cartesian) {
        float motors[MAX_N_AXIS];
        
        if (!transform_cartesian_to_motors(motors, cartesian)) {
            log_info("SCARA soft limit at (" << cartesian[X_AXIS] << ", " << cartesian[Y_AXIS] << ")");
            limit_error();
            return true;
        }
        
        return false;
    }

    /*
      Constrain jog moves to valid work envelope
    */
    void Scara::constrain_jog(float* target, plan_line_data_t* pl_data, float* position) {
        float motor_pos[MAX_N_AXIS] = { 0.0f };
        
        // If the target is reachable, do nothing
        if (transform_cartesian_to_motors(motor_pos, target)) {
            return;
        }
        
        log_warn("SCARA kinematics soft limit - jog rejected");
        copyAxes(target, position);
        
        // TODO: Could walk back from target in small increments to find a valid position
    }

    /*
      Check if homing is possible
    */
    bool Scara::canHome(AxisMask axisMask) {
        if (ambiguousLimit()) {
            log_error("Ambiguous limit switch touching. Manually clear all switches");
            return false;
        }
        return true;
    }

    /*
      Release motors that are still needed for homing.
      For SCARA, each axis/motor is independent in homing (joint space),
      so we just unlimit each motor individually.
    */
    void Scara::releaseMotors(AxisMask axisMask, MotorMask motors) {
        auto n_axis = Axes::_numberAxis;
        for (axis_t axis = X_AXIS; axis < n_axis; axis++) {
            if (bitnum_is_true(axisMask, axis)) {
                Stepping::unlimit(axis, MOTOR0);
            }
        }
    }

    /*
      Handle limit switch events during homing.
      Each SCARA joint homes independently, so when a limit is hit
      we simply clear that motor from the mask.
    */
    bool Scara::limitReached(AxisMask& axisMask, MotorMask& motors, MotorMask limited) {
        // Clear the motors whose limits have been reached
        clear_bits(motors, limited);

        auto oldAxisMask = axisMask;

        // Set axisMask according to the motors that are still running
        axisMask = Machine::Axes::motors_to_axes(motors);

        // Return true when an axis drops out of the mask, causing replan
        return axisMask != oldAxisMask;
    }

    /*
      SCARA homing move - operates directly in motor/joint space.
      
      Unlike normal cartesian_to_motors() which does inverse kinematics,
      homing sends motor positions (angles in degrees) directly via
      mc_move_motors(), bypassing the kinematic transform.
      
      This ensures that $HX only moves the theta motor and
      $HY only moves the psi motor.
      
      The implementation reuses axesVector() from Cartesian to compute
      distance and speed, but uses get_motor_pos() as the starting point
      (motor angles in degrees) instead of get_mpos() (cartesian coords).
      The resulting target is sent directly to mc_move_motors() without
      going through inverse kinematics.
    */
    void Scara::homing_move(AxisMask axisMask, MotorMask motors, Machine::Homing::Phase phase, uint32_t settling_ms) {
        releaseMotors(axisMask, motors);

        auto axes   = config->_axes;
        auto n_axis = axes->_numberAxis;

        // Use axesVector logic but in motor (angle) space instead of cartesian space.
        // axesVector computes: target = startPos + distance, rate.
        // We compute the same thing but with motor positions as the starting point.
        
        float rate;
        float target[MAX_N_AXIS];
        
        // Cartesian::axesVector uses get_mpos() internally as start position.
        // We need motor positions instead, so we call axesVector and then
        // replace the result with motor-space calculations.
        
        // Start from current motor positions (angles in degrees)
        copyAxes(target, get_motor_pos());
        
        float maxSeekTime = 0.0f;
        float ratesq      = 0.0f;
        settling_ms       = 0;

        float rates[MAX_N_AXIS]    = { 0 };
        float distance[MAX_N_AXIS] = { 0 };

        bool seeking  = phase == Machine::Homing::Phase::FastApproach;
        bool approach = seeking || phase == Machine::Homing::Phase::SlowApproach;

        for (axis_t axis = X_AXIS; axis < n_axis; axis++) {
            if (bitnum_is_false(axisMask, axis)) {
                continue;
            }

            auto axisConfig = axes->_axis[axis];
            auto homing     = axisConfig->_homing;
            if (!homing) {
                continue;
            }

            settling_ms = std::max(settling_ms, homing->_settle_ms);

            float axis_rate = 1;
            float travel    = 0;
            switch (phase) {
                case Machine::Homing::Phase::FastApproach:
                    axis_rate = homing->_seekRate;
                    travel    = axisConfig->_maxTravel;
                    break;
                case Machine::Homing::Phase::PrePulloff:
                case Machine::Homing::Phase::SlowApproach:
                case Machine::Homing::Phase::Pulloff0:
                case Machine::Homing::Phase::Pulloff1:
                    axis_rate = homing->_feedRate;
                    travel    = axisConfig->commonPulloff();
                    break;
                case Machine::Homing::Phase::Pulloff2:
                    axis_rate = homing->_feedRate;
                    travel    = axisConfig->extraPulloff();
                    break;
                default:
                    break;
            }

            switch (phase) {
                case Machine::Homing::Phase::PrePulloff: {
                    MotorMask axisMotors = Machine::Axes::axes_to_motors(1 << axis);
                    bool      posLimited = bits_are_true(Machine::Axes::posLimitMask, axisMotors);
                    bool      negLimited = bits_are_true(Machine::Axes::negLimitMask, axisMotors);
                    if (posLimited) {
                        distance[axis] = -travel;
                    } else if (negLimited) {
                        distance[axis] = travel;
                    } else {
                        distance[axis] = 0;
                    }
                } break;

                case Machine::Homing::Phase::FastApproach:
                case Machine::Homing::Phase::SlowApproach:
                    distance[axis] = homing->_positiveDirection ? travel : -travel;
                    break;

                case Machine::Homing::Phase::Pulloff0:
                case Machine::Homing::Phase::Pulloff1:
                case Machine::Homing::Phase::Pulloff2:
                    distance[axis] = homing->_positiveDirection ? -travel : travel;
                    break;

                default:
                    break;
            }

            ratesq += (axis_rate * axis_rate);
            rates[axis] = axis_rate;

            auto seekTime = travel / axis_rate;
            if (seekTime > maxSeekTime) {
                maxSeekTime = seekTime;
            }
        }

        // Scale distances and apply to target
        for (axis_t axis = X_AXIS; axis < n_axis; axis++) {
            if (bitnum_is_false(axisMask, axis)) {
                continue;
            }
            if (phase == Machine::Homing::Phase::FastApproach) {
                float absDistance = maxSeekTime * rates[axis];
                distance[axis]   = distance[axis] >= 0 ? absDistance : -absDistance;
            }

            auto axisConfig = axes->_axis[axis];
            auto homing     = axisConfig->_homing;
            if (homing) {
                auto scaler = approach ? (seeking ? homing->_seek_scaler : homing->_feed_scaler) : 1.0;
                distance[axis] *= scaler;
                target[axis] += distance[axis];
            }
        }

        rate = sqrtf(ratesq);
        if (rate == 0) {
            return;
        }

        log_debug("SCARA homing motor target " << target[0] << "," << target[1] << "," << target[2] << " @ " << rate);

        plan_line_data_t plan_data      = {};
        plan_data.spindle_speed         = 0;
        plan_data.motion                = {};
        plan_data.motion.systemMotion   = 1;
        plan_data.motion.noFeedOverride = 1;
        plan_data.spindle               = SpindleState::Disable;
        plan_data.coolant               = {};
        plan_data.line_number           = 0;
        plan_data.is_jog                = false;
        plan_data.feed_rate             = rate;

        mc_move_motors(target, &plan_data);
        protocol_send_event(&cycleStartEvent);
    }

    /*
      Set the machine position after homing.
      
      Called by Homing::set_mpos() after limit switches are triggered.
      Homing::set_mpos() sets mpos[axis] = homing->_mpos for each homed axis.
      
      For SCARA X/Y axes, homing operates in joint/motor space (degrees).
      The homing->_mpos value is the motor angle (degrees) at the limit switch.
      We set the motor steps directly from these angle values WITHOUT
      doing inverse kinematics - because homing already knows the motor
      angles, not cartesian positions.
      
      For other axes (Z servo etc.), the values pass through unchanged.
      
      After setting motor steps, we update mpos[] to the corresponding
      cartesian coordinates via forward kinematics, so the rest of the
      system (GCode parser, planner) sees correct cartesian positions.
    */
    void Scara::set_homed_mpos(float* mpos) {
        auto  n_axis = Axes::_numberAxis;

        if (_calibrating) {
            // CALIBRATION MODE:
            // At this point, the motor has physically moved from the calibration
            // pose to (limit switch - pulloff). The step counters reflect the
            // actual position because we set them to the init angles before homing.
            //
            // Read the current steps and convert to angles - this gives us
            // the actual physical angle after homing (at the pulloff position).
            //
            // However, mpos_mm in the homing config represents the position
            // that the system assigns after homing completes (after pulloff).
            // So we save the current angle as the mpos_mm value.
            
            float* current_motor_pos = get_motor_pos();
            _cal_result_theta = current_motor_pos[X_AXIS];
            _cal_result_psi   = current_motor_pos[Y_AXIS];

            log_debug("Calibration captured angles: theta=" << _cal_result_theta 
                      << " psi=" << _cal_result_psi);

            // Still set the steps correctly for the current position
            // (using the angles we just read, which are already in the step counter)
            // and update cartesian position via forward kinematics
            float motor_pos[MAX_N_AXIS];
            for (axis_t axis = X_AXIS; axis < n_axis; axis++) {
                motor_pos[axis] = current_motor_pos[axis];
            }
            motors_to_cartesian(mpos, motor_pos, n_axis);
            return;
        }

        // NORMAL HOMING MODE:
        // For SCARA, mpos[X] and mpos[Y] from Homing::set_mpos() are
        // motor angles (degrees) at the limit switch position.
        // Set motor steps directly from these angle values.
        float motor_pos[MAX_N_AXIS];
        for (axis_t axis = X_AXIS; axis < n_axis; axis++) {
            motor_pos[axis] = mpos[axis];
        }

        for (axis_t axis = X_AXIS; axis < n_axis; axis++) {
            set_steps(axis, motor_pos_to_steps(motor_pos[axis], axis));
        }

        // Now update mpos[] to cartesian coordinates via forward kinematics
        // so that the GCode parser and planner see correct positions.
        // motors_to_cartesian converts motor angles → cartesian for X/Y
        // and passes through Z and other axes unchanged.
        motors_to_cartesian(mpos, motor_pos, n_axis);
    }

    /*
      Auto-calibration procedure ($M700)
      
      Prerequisites:
        User must manually position the arm at the calibration pose:
        - theta (upper arm) = _cal_init_theta (default 0°, along -X axis)
        - psi (lower arm)   = _cal_init_psi   (default 45°)
      
      Algorithm:
        1. Set step counters to correspond to the calibration pose angles
        2. Enable calibration mode flag
        3. Run standard homing - arms rotate to limit switches, then pulloff
        4. In set_homed_mpos (calibration mode), capture the current motor 
           angles from step counters BEFORE they get overwritten
        5. These captured angles become the new mpos_mm values
        6. Update homing config and re-sync position
        
      After calibration, user should save with $Config/Save (NVS).
      Subsequent power-ups only need G28 ($H) to home.
    */
    Error Scara::auto_calibrate(Channel& out) {
        if (!state_is(State::Idle)) {
            log_error("Cannot calibrate: machine not idle");
            return Error::IdleError;
        }

        auto axes   = config->_axes;
        auto n_axis = axes->_numberAxis;

        log_info("SCARA auto-calibration starting...");
        log_info("  Initial pose: theta=" << _cal_init_theta << " psi=" << _cal_init_psi);

        // Step 1: Set motor position to the known calibration pose angles.
        // This establishes the step counter baseline so that when homing moves
        // the motors, the step counters track the actual physical angles.
        set_motor_pos(X_AXIS, _cal_init_theta);
        set_motor_pos(Y_AXIS, _cal_init_psi);

        // Update cartesian position via forward kinematics for system consistency
        float motor_pos[MAX_N_AXIS];
        motor_pos[X_AXIS] = _cal_init_theta;
        motor_pos[Y_AXIS] = _cal_init_psi;
        for (axis_t axis = Z_AXIS; axis < n_axis; axis++) {
            motor_pos[axis] = steps_to_motor_pos(get_axis_steps(axis), axis);
        }
        float cartesian[MAX_N_AXIS];
        motors_to_cartesian(cartesian, motor_pos, n_axis);
        gc_sync_position();

        log_info("  Motor position set to calibration pose");

        // Step 2: Enable calibration mode. In this mode, set_homed_mpos() will
        // capture the current motor angles (from step counters) instead of
        // overwriting them with the old mpos_mm values.
        _calibrating = true;
        _cal_result_theta = 0;
        _cal_result_psi = 0;

        log_info("  Running homing sequence...");

        // Step 3: Run homing for X and Y axes (cycle 1).
        // Homing moves motors to limit switches, does pulloff, then calls
        // set_homed_mpos(). In calibration mode, set_homed_mpos captures
        // the angle at the post-pulloff position.
        AxisMask xy_mask = bitnum_to_mask(X_AXIS) | bitnum_to_mask(Y_AXIS);
        Machine::Homing::run_cycles(xy_mask);

        // Wait for homing to complete
        do {
            protocol_execute_realtime();
        } while (state_is(State::Homing));

        _calibrating = false;

        if (state_is(State::Alarm)) {
            log_error("Calibration failed: homing resulted in alarm");
            return Error::InvalidStatement;
        }

        log_info("  Homing complete");
        log_info("  Calibrated angles: theta=" << _cal_result_theta << " psi=" << _cal_result_psi);

        // Step 4: Update the homing mpos_mm config with the calibrated values.
        // These values represent the motor angles at the post-pulloff position
        // after homing. On future G28 commands, set_homed_mpos will use these
        // values to correctly set the motor position.
        auto x_homing = axes->_axis[X_AXIS]->_homing;
        auto y_homing = axes->_axis[Y_AXIS]->_homing;

        if (x_homing) {
            x_homing->_mpos = _cal_result_theta;
            log_info("  X axis homing mpos_mm = " << _cal_result_theta);
        }
        if (y_homing) {
            y_homing->_mpos = _cal_result_psi;
            log_info("  Y axis homing mpos_mm = " << _cal_result_psi);
        }

        // Step 5: Re-run set_homed_mpos in normal mode with the correct values
        // to ensure the step counters and cartesian position are fully consistent.
        float mpos[MAX_N_AXIS];
        float* current_mpos = get_mpos();
        for (axis_t axis = X_AXIS; axis < n_axis; axis++) {
            mpos[axis] = current_mpos[axis];
        }
        mpos[X_AXIS] = _cal_result_theta;
        mpos[Y_AXIS] = _cal_result_psi;
        set_homed_mpos(mpos);

        gc_sync_position();

        // Step 6: Save the calibration init pose (e.g. theta=0, psi=45) as G28 position.
        // After homing ($H), the user can send G28 to move back to this known pose.
        // Calculate the cartesian coordinates of the init pose via forward kinematics.
        {
            float init_motor[MAX_N_AXIS] = {};
            init_motor[X_AXIS] = _cal_init_theta;
            init_motor[Y_AXIS] = _cal_init_psi;
            for (axis_t axis = Z_AXIS; axis < n_axis; axis++) {
                init_motor[axis] = 0;
            }
            float init_cartesian[MAX_N_AXIS] = {};
            motors_to_cartesian(init_cartesian, init_motor, n_axis);

            coords[CoordIndex::G28]->set(init_cartesian);
            gc_ngc_changed(CoordIndex::G28);

            log_info("  G28 position saved: X=" << init_cartesian[X_AXIS] 
                     << " Y=" << init_cartesian[Y_AXIS]);
        }

        log_info("SCARA calibration complete!");
        log_msg_to(out, "Calibration complete. Limit angles: theta=" 
                   << _cal_result_theta << " psi=" << _cal_result_psi);
        log_msg_to(out, "Run $CD=/localfs/config.yaml to save config.");

        return Error::Ok;
    }

    // Configuration registration
    namespace {
        KinematicsFactory::InstanceBuilder<Scara> registration("Scara");
    }

}  // namespace Kinematics
