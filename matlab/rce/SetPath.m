% Script to set up path to use the tc_seed routines

restoredefaultpath;
addpath('~/etc/matlab/common', '~/etc/matlab/m_map', '~/etc/matlab/metpack', '~/etc/matlab/obj_anal', '~/etc/matlab/rce', '-end');
run ~/etc/matlab/nctoolbox/setup_nctoolbox.m;