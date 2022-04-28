
;;***** Inputs
path_bkg_file                            = '{{ background_L1A_fits_filename }}'
path_sci_file                            = '{{ flaring_data_L1A_fits_filename }}'
flare_start_UTC                          = '{{ flare_start_UTC }}'
flare_end_UTC                            = '{{ flare_end_UTC }}'
energy_range_science_channel_upper_limit = '{{ energy_range_science_channel_upper_limit }}'
energy_range_science_channel_lower_limit = '{{ energy_range_science_channel_lower_limit }}'
vis_fwdfit_map_filename                  = '{{ vis_fwdfit_map_filename }}'
bp_map_filename                          = '{{ bp_map_filename }}'
L0                                       = '{{ L0 }}'
B0                                       = '{{ B0 }}'
RSUN                                     = '{{ apparent_radius_sun }}'
roll_angle                               = '{{ roll_angle_solo }}'

;;***** Parameters
energy_range = [energy_range_science_channel_lower_limit,energy_range_science_channel_upper_limit]
mapcenter    = [0., 0.]
subc_index   = stix_label2ind(['10a','10b','10c','9a','9b','9c','8a','8b','8c','7a','7b','7c'])
silent       = 1
imsize_bp    = [512,512]
pixel_bp     = RSUN * 2.6 / imsize_bp
imsize       = [256,256]
pixel        = [1.,1.]
time_range   = anytim([flare_start_UTC,flare_end_UTC])

;;******* Compute the visibility values
vis  = stix2vis_sep2021(path_sci_file, time_range, energy_range, mapcenter, subc_index=subc_index, $
                        path_bkg_file=path_bkg_file, silent=silent)

;;******* Use Giordano's (u,v) points: no need to perform projection correction (see Giordano et al., 2015)
subc_str = stx_construct_subcollimator()
uv       = stx_uv_points_giordano()
u        = -uv.u * subc_str.phase
v        = -uv.v * subc_str.phase
vis.u    = u[subc_index]
vis.v    = v[subc_index]

;;******* Compute the Back Projection map
bp_nat_map = stx_bproj(vis,imsize_bp,pixel_bp)

;;******* Compute the coordinates of the maximum value of the Back Projection map, i.e. of the location of the flare
max_bp       = max(bp_nat_map.data, ind_max)
ind_max      = array_indices(bp_nat_map.data, ind_max)
max_bp_coord = [(ind_max[0]-imsize_bp[0]/2)*pixel_bp[0]+mapcenter[0], (ind_max[1]-imsize_bp[1]/2)*pixel_bp[1]+mapcenter[1]]

;;***** Re-compute the visibilities (needed for setting the map center correctly and for using subcollimators 3 to 10)
vis  = stix2vis_sep2021(path_sci_file, time_range, energy_range, max_bp_coord, $
                        xy_flare=max_bp_coord, path_bkg_file=path_bkg_file, silent=silent)

;;******* Compute the FWDFIT solution
vis_fwdfit_pso_map = stx_vis_fwdfit_pso('circle', vis, param_opt=param_opt, imsize=imsize, pixel=pixel, $
                                        srcstr = srcstrout_pso, silent=silent)
vis_fwdfit_pso_map.L0 = L0
vis_fwdfit_pso_map.B0 = B0
vis_fwdfit_pso_map.RSUN = RSUN
vis_fwdfit_pso_map.roll_angle = roll_angle

map2fits, bp_nat_map, bp_map_filename
map2fits, vis_fwdfit_pso_map, vis_fwdfit_map_filename
