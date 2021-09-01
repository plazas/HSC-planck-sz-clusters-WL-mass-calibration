import os
import glob
import yaml

root="/project/plazas/WORK/HSC/weaklens_pipeline/"
list_randoms = glob.glob(root+"DataStore/randoms_S19A/*")
# /project/plazas/WORK/HSC/weaklens_pipeline/DataStore/randoms_S19A/random_points_s19a_98.dat
config_template = root + "configs/config.andres-random-centers-p-cuts-2021JUN04"
config_dir = root + "configs/configs_100randoms_boost_factor/"
wl_code = root + "weaklens_aroundsource.py"

for randoms_file in list_randoms[2:]:
    print (f"Processing file {randoms_file}")
    with open(config_template, 'r') as yamlfile:
        config = yaml.safe_load(yamlfile)
    
    out_dir_random = randoms_file.split('/')[-1].split('.')[0]
    out_dir = root + "output_WL/100RANDOMS_2021JUL01/" + out_dir_random

    # replace entries that change
    config['outputdir'] = out_dir
    config['lens']['fname'] = randoms_file

    # save new config in new file
    new_config_file_name = config_dir + f"config.{out_dir_random}-p-cuts"
    with open(new_config_file_name, 'w') as yamlfile:
        yaml.safe_dump(config, yamlfile)
    
    # Run the weak lensing code
    cmd = f"python {wl_code} --config {new_config_file_name}"
    print ("Running: ")
    print (cmd)
    os.system(cmd)
