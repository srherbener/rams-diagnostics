% Script to set up path to use the ideal_etc routines

restoredefaultpath;
addpath('~/etc/matlab/common', '~/etc/matlab/m_map', '~/etc/matlab/obj_anal', '~/etc/matlab/ideal_etc', '~/etc/matlab/metpack', '-end');

run ~/etc/matlab/nctoolbox/setup_nctoolbox.m;
