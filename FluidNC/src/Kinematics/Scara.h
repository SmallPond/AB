// Copyright (c) 2026 -  dingmos
// Use of this source code is governed by a GPLv3 license that can be found in the LICENSE file.

#pragma once

/*
	Scara.h

	This implements SCARA (Selective Compliance Assembly Robot Arm) Kinematics
	with absolute angle mode.

	In absolute angle mode:
	- The first arm angle (theta) is measured from the -X axis
	- The second arm angle (psi) is measured from the world coordinate system,
	  not relative to the first arm

	This is different from "relative angle" SCARA where the second arm angle
	is measured relative to the first arm.

	References:
	- https://en.wikipedia.org/wiki/SCARA
*/

#include "Kinematics.h"
#include "Cartesian.h"

// M_PI is not defined in standard C/C++ but some compilers
// support it anyway.  The following suppresses Intellisense
// problem reports.
#ifndef M_PI
#    define M_PI 3.14159265358979323846
#endif

namespace Kinematics {

    class Scara : public Cartesian {
    public:
        Scara(const char* name) : Cartesian(name) {}

        Scara(const Scara&)            = delete;
        Scara(Scara&&)                 = delete;
        Scara& operator=(const Scara&) = delete;
        Scara& operator=(Scara&&)      = delete;

        // Kinematic Interface
        virtual void init() override;
        virtual void init_position() override;
        
        bool cartesian_to_motors(float* target, plan_line_data_t* pl_data, float* position) override;
        void motors_to_cartesian(float* cartesian, float* motors, axis_t n_axis) override;
        bool transform_cartesian_to_motors(float* motors, float* cartesian) override;
        
        virtual bool invalid_line(float* cartesian) override;
        virtual void constrain_jog(float* cartesian, plan_line_data_t* pl_data, float* position) override;
        
        bool canHome(AxisMask axisMask) override;
        void releaseMotors(AxisMask axisMask, MotorMask motors) override;
        bool limitReached(AxisMask& axisMask, MotorMask& motors, MotorMask limited) override;
        void homing_move(AxisMask axes, MotorMask motors, Machine::Homing::Phase phase, uint32_t settling_ms) override;
        void set_homed_mpos(float* mpos) override;

        // Configuration handlers
        void validate() override {}
        virtual void group(Configuration::HandlerBase& handler) override;
        void afterParse() override {}

        // Calibration support
        // Called by $M700 command to perform auto-calibration
        // The arm must be pre-positioned at the calibration pose:
        //   theta = _cal_init_theta (default 0°), psi = _cal_init_psi (default 45°)
        // Steps: 1) Record current step counts
        //        2) Run homing to hit limit switches
        //        3) Calculate limit switch angles from initial pose + moved steps
        //        4) Update homing mpos_mm with calculated angles
        Error auto_calibrate(Channel& out) override;

        ~Scara() {}

    private:
        // Forward kinematics: Convert motor angles (in degrees) to cartesian coordinates
        // Angles are in degrees, cartesian is in mm
        void forward_kinematics(float* cartesian, float theta_deg, float psi_deg);
        
        // Inverse kinematics: Convert cartesian coordinates to motor angles
        // Returns true if the position is reachable, false otherwise
        bool inverse_kinematics(float cartesian_x, float cartesian_y, float& theta_deg, float& psi_deg);

        // Helper functions for angle conversions
        inline float degrees_to_radians(float degrees) { return degrees * (M_PI / 180.0f); }
        inline float radians_to_degrees(float radians) { return radians * (180.0f / M_PI); }

        // Configuration parameters
        float _linkage1_mm = 100.0f;  // Length of first (upper) arm in mm
        float _linkage2_mm = 100.0f;  // Length of second (lower) arm in mm
        
        // Offset from world coordinate system to user coordinate system
        float _scara_offset_x = 0.0f;
        float _scara_offset_y = 0.0f;
        
        // Segment length for breaking up non-linear moves
        float _kinematic_segment_len_mm = 1.0f;

        // Calibration initial pose (degrees)
        // Default: upper arm along -X axis (0°), elbow at 45°
        float _cal_init_theta = 0.0f;
        float _cal_init_psi   = 45.0f;

        // Calibration runtime state
        bool  _calibrating      = false;
        float _cal_result_theta = 0.0f;
        float _cal_result_psi   = 0.0f;

        // Motor position storage for feedrate calculation
        float _last_motor_pos[MAX_N_AXIS] = { 0 };
    };

}  // namespace Kinematics
