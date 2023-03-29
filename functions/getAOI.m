function [chansOTleft,chansOTright,chansOTall] = getAOI()


%% MEG Area of Interests 

 sensors = neuromagGetSensorAOIs;

 % Left AOIs
 megAOI_left = unique([sensors.frontal.left, sensors.parietal.left,...
                     sensors.temporal.left, sensors.occipital.left])';
 megAOI_left_label = meg.header.label(megAOI_left);

 % Right AOIs
 megAOI_right = unique([sensors.frontal.right, sensors.parietal.right,...
                      sensors.temporal.right, sensors.occipital.right])';
 megAOI_right_label = meg.header.label(megAOI_right);

 % Left and Right AOIs
 megAOI_All = unique([sensors.frontal.left, sensors.parietal.left,...
                    sensors.temporal.left, sensors.occipital.left,...
                    sensors.frontal.right, sensors.parietal.right,...
                    sensors.temporal.right, sensors.occipital.right])';
 megAOI_all_label = meg.header.label(megAOI_All);

 %% channel of interest (occipital)

 chansOleft = unique([sensors.occipital.left])';
 chansOleft(chansOleft==80) = [];

 chansOright = unique([sensors.occipital.right])';
 chansOright(chansOright==79) = [];

 chansOall = unique([sensors.occipital.left, sensors.occipital.right])';

 %% channel of interest (temporal)

 chansTleft = unique([sensors.temporal.left])';
 chansTright = unique([sensors.temporal.right])';
 chansTall = unique([sensors.temporal.left, sensors.temporal.right])';

%% channel of interest (occipito-temporal)

chansOTleft = unique(cat(1,chansOleft,chansTleft));
chansOTright = unique(cat(1,chansOright,chansTright));
chansOTall = unique(cat(1,chansOall,chansTall));

end