
folder='/Users/admin/Documents/GitHub/stix_quicklook_images/'
data_folder=folder+'/data/'
if ~file_exist(data_folder) then file_mkdir, data_folder
fits_files_folder=folder+'/fits files/'
if ~file_exist(fits_files_folder) then file_mkdir, fits_files_folder

;; Download data
website_url = 'https://pub023.cs.technik.fhnw.ch/download/fits/bsd/'
uid = '2105070034'
sock_copy, website_url + uid, out_name, status = status, out_dir = data_folder, local_file=path_sci_file, clobber=0
uid = '2105080012'
sock_copy, website_url + uid, out_name, status = status, out_dir = data_folder, local_file=path_bkg_file, clobber=0


flare_start_UTC='7-May-2021 18:51:00'
flare_end_UTC  ='7-May-2021 18:53:40'
energy_range_full_disk_bp_map_lower_limit_keV=6
energy_range_full_disk_bp_map_upper_limit_keV=10
energy_range_lower_limit_keV=22
energy_range_upper_limit_keV=50
full_disk_bp_map_filename=fits_files_folder + 'full_disk_bp_map.fits' 
full_disk_bp_map_size=[512,512]
full_disk_bp_map_subc_index=stix_label2ind(['10a','10b','10c','9a','9b','9c','8a','8b','8c','7a','7b','7c'])
full_disk_bp_map_mapcenter=[0.,0.]
map_size=[256,256]
pixel_size=[1.,1.]
subc_index=stix_label2ind(['10a','10b','10c','9a','9b','9c','8a','8b','8c','7a','7b','7c',$
                          '6a','6b','6c','5a','5b','5c','4a','4b','4c','3a','3b','3c'])
bp_map_filename=fits_files_folder + 'bp_map.fits' 
vis_fwdfit_map_filename=fits_files_folder + 'vis_fwdfit_map.fits' 
vis_fwdfit_source_type='multi'
em_map_filename=fits_files_folder + 'em_map.fits'
clean_map_filename=fits_files_folder + 'clean_map.fits'
clean_niter  = 200    
clean_gain   = 0.1    
clean_beam_width = 20.
clean_uniform_weighting=0
L0=-97.444793247335838
B0=-0.31161748599820238
RSUN=1045.9928616497800
roll_angle=-1.9719335530048476
x_offset_arcsec=0
y_offset_arcsec=0    


stx_image_reconstruct, path_bkg_file, path_sci_file, $
                        flare_start_UTC, flare_end_UTC, $
                        energy_range_lower_limit_keV, energy_range_upper_limit_keV, $
                        energy_range_full_disk_bp_map_lower_limit_keV, energy_range_full_disk_bp_map_upper_limit_keV, $
                        full_disk_bp_map_filename, full_disk_bp_map_size, full_disk_bp_map_subc_index, full_disk_bp_map_mapcenter, $
                        map_size, pixel_size, subc_index, $
                        bp_map_filename, $
                        vis_fwdfit_map_filename, vis_fwdfit_source_type, $
                        em_map_filename, $
                        clean_map_filename, clean_niter, clean_gain, clean_beam_width, clean_uniform_weighting, $
                        L0, B0, RSUN, roll_angle, $
                        x_offset_arcsec, y_offset_arcsec     

fits2map,full_disk_bp_map_filename,full_disk_bp_map
fits2map,bp_map_filename,bp_map
fits2map,em_map_filename,em_map
fits2map,vis_fwdfit_map_filename,vis_fwdfit_map
fits2map,clean_map_filename, clean_map

chs2=1.
window,0,xsize=500,ysize=500
cleanplot
loadct, 3
plot_map,full_disk_bp_map,charsize=chs2,title='Full disk Back Projection', /limb,grid_spacing=15

loadct, 5
window,1,xsize=800,ysize=800
!p.multi=[0,2,2]
plot_map,clean_map,charsize=chs2,title='CLEAN', /limb,grid_spacing=15
plot_map,bp_map,charsize=chs2,title='Back Projection',/limb,grid_spacing=5
plot_map,em_map,charsize=chs2,title='EM',/limb,grid_spacing=5
plot_map,vis_fwdfit_map,charsize=chs2,title='VIS_FWDFIT_PSO',/limb,grid_spacing=5


end