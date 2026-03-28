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
        // SCARA can home using standard cartesian-style homing
        // The motors are directly controlled for homing
        return Cartesian::canHome(axisMask);
    }

    // Configuration registration
    namespace {
        KinematicsFactory::InstanceBuilder<Scara> registration("Scara");
    }

}  // namespace Kinematics
