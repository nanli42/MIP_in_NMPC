#!/bin/bash

method_opt=1 # 1 - MIP
track_opt=$1 # 1 - track1, 2 - track2

single_RC_mode_opt=0 # 0 - no, 1 - yes (for testing single-vehicle racing mode)
LV_traj_generation_opt=0 # 0 - no, 1 - yes (for generating LV's trajectory)

# change track
cmd_string="%s/#define TRACK_OPT .*/#define TRACK_OPT ${track_opt}/g | wq";
vim -c "$cmd_string" ../src/track.cpp
cmd_string="%s/#define TRACK_OPT .*/#define TRACK_OPT ${track_opt}/g | wq";
vim -c "$cmd_string" ../src/plot.cpp
# cmd_string="%s/#define TRACK_OPT .*/#define TRACK_OPT ${track_opt}/g | wq";
# vim -c "$cmd_string" ../src/shadow_vehicle.cpp
cmd_string="%s/#define TRACK_OPT .*/#define TRACK_OPT ${track_opt}/g | wq";
vim -c "$cmd_string" ../include/LV_traj.h

# change definition in main.cpp
cmd_string="%s/#define GENERATE_LV_TRAJ .*/#define GENERATE_LV_TRAJ ${LV_traj_generation_opt}/g |"
cmd_string+="%s/#define singRC .*/#define singRC ${single_RC_mode_opt}/g |"
cmd_string+="%s/#define METHOD .*/#define METHOD ${method_opt}/g |"
cmd_string+="wq";
vim -c "$cmd_string" ../main.cpp

make
