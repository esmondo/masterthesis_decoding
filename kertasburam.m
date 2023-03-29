close all
clear all
clc

addpath ~/any_results

dat2mit = load('dat2_navg1_megmitrej.mat');
dat2ohne = load('dat2_navg1_megohnerej.mat');

dat4ohne = load('dat4_navg1_megeegohnerej.mat');

Cmit = cat(1, dat2mit.akurasi_meg, dat2mit.akurasiSTF_meg);
Cohne = cat(1, dat2ohne.akurasi_meg, dat2ohne.akurasiSTF_meg);

% Cdiff = Cmit - Cohne;