% Script to set up path to use the tc_seed routines

restoredefaultpath;
addpath('~/etc/matlab/common', '~/etc/matlab/m_map', '~/etc/matlab/metpack', '~/etc/matlab/obj_anal', '~/etc/matlab/pivmat', '~/etc/matlab/ts_debby', '-end');
run ~/etc/matlab/nctoolbox/setup_nctoolbox.m;
