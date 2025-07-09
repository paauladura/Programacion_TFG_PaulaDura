%-------------------------------------------------------------------------------
% Add path for SHbundle, gracebundle, uberall software packages 
addpath(genpath('/Users/paula/OneDrive/Escritorio/4º CUATRI 2/TFG/CosasApasar'))
%
%-------------------------------------------------------------------------------
% 
%% Import the static gravity field model ITSG-Grace2018s with readicgem. 
   [scs, ncs, header, scst, ncst]  = readicgem('ITSG-Grace2018s.gfc');
  
%-------------------------------------------------------------------------------

%% Geoid Heights
%
% Compute global grid of geoid heights up to various maximum harmonic degrees 
% with gshs_. 
 [geoid_f_10, theRAD_10, lamRAD_10] = gshs_(dis,'quant','geoid','max_lm',10,'sub_WGS84',false);
 [geoid_f_50, theRAD_50, lamRAD_50] = gshs_(dis,'quant','geoid','max_lm',50,'sub_WGS84',false);
 [geoid_f_100, theRAD_100, lamRAD_100] = gshs_(dis,'quant','geoid','max_lm',100,'sub_WGS84',false);
 [geoid_f_200, theRAD_200, lamRAD_200] = gshs_(dis,'quant','geoid','max_lm',200,'sub_WGS84',false);
 
%% Plot the results 
mapfield(geoid_f_10,'pole',1);
title('n=10', 'FontSize',30)
colormap('jet')
cb=colorbar()
cb.Label.String = 'Geoide';
cb.Label.FontSize = 27;
cb.Label.FontWeight = 'bold';
set(gca, 'FontSize', 20) ;
 
 mapfield(geoid_f_50,'pole',1);
title('n=50', 'FontSize',30)
colormap(jet)
cb=colorbar()
cb.Label.String = 'Geoide';
cb.Label.FontSize = 27;
cb.Label.FontWeight = 'bold';
set(gca, 'FontSize', 20) ;
 
mapfield(geoid_f_100,'pole',1);
title('n=100', 'FontSize',30)
colormap(jet)
cb=colorbar()
cb.Label.String = 'Geoide';
cb.Label.FontSize = 27;
cb.Label.FontWeight = 'bold';
set(gca, 'FontSize', 20) ;
 
mapfield(geoid_f_200,'pole',1);
title('n=200', 'FontSize',30)
colormap(jet)
cb=colorbar()
cb.Label.String = 'Geoide';
cb.Label.FontSize = 27;
cb.Label.FontWeight = 'bold';
set(gca, 'FontSize', 20) ;
 % As the degree of the harmonics go up, the results are more detailed.
 

%% Potential
%
% Compute the potential up to nmax = 200 as seen from different 
% satellite altitudes (e.g. h = 250, 450, 1000, 20000km). 
[tr_f_200_h_250, theRAD_200, lamRAD_200] = gshs_(dis,'quant','potential','max_lm',200,...
                                      'sub_WGS84',false,'height',250*1000);

[tr_f_200_h_450, theRAD_200, lamRAD_200] = gshs_(dis,'quant','potential','max_lm',200,...
                                      'sub_WGS84',false,'height',450*1000);
                                  
[tr_f_200_h_10k, theRAD_200, lamRAD_200] = gshs_(dis,'quant','potential','max_lm',200,...
                                      'sub_WGS84',false,'height',1000*1000);
                                   
 [tr_f_200_h_20k theRAD_200, lamRAD_200] = gshs_(dis,'quant','potential','max_lm',200,...
                                      'sub_WGS84',false,'height',250*20000);  
                                  
%% Visualize the results in terms of global grids 
mapfield(tr_f_200_h_250,'pole',1);
title('h = 250km', 'FontSize',30)
colormap('jet')
cb=colorbar()
cb.Label.String = 'Potencial (V)';
cb.Label.FontSize = 27;
cb.Label.FontWeight = 'bold';
set(gca, 'FontSize', 20) ;

mapfield(tr_f_200_h_450,'pole',1);
title('h = 450km' ,'FontSize',30)
colormap('jet')
cb=colorbar()
cb.Label.String = 'Potencial (V)';
cb.Label.FontSize = 27;
cb.Label.FontWeight = 'bold';
set(gca, 'FontSize', 20) ;

mapfield(tr_f_200_h_10k,'pole',1);
title('h = 10.000km' ,'FontSize',30)
colormap('jet')
cb=colorbar()
cb.Label.String = 'Potencial (V)';
cb.Label.FontSize = 27;
cb.Label.FontWeight = 'bold';
set(gca, 'FontSize', 20) ;

mapfield(tr_f_200_h_20k,'pole',1);
title('h = 20.000km' ,'FontSize',30)
colormap('jet')
cb=colorbar()
cb.Label.String = 'Potencial (V)';
cb.Label.FontSize = 27;
cb.Label.FontWeight = 'bold';
set(gca, 'FontSize', 20) ;

%-------------------------------------------------------------------------------
