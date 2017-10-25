% Script to set up path to use the tempest routines

restoredefaultpath;
addpath('~/etc/matlab/common', '~/etc/matlab/m_map', '~/etc/matlab/obj_anal', '~/etc/matlab/tempest', '~/etc/matlab/metpack', '-end');

run ~/etc/matlab/nctoolbox/setup_nctoolbox.m;
startup
