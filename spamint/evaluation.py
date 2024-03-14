import subprocess
import os


def runSpaTalk(rscript_executable = '/apps/software/R/4.2.0-foss-2021b/bin/Rscript',
               st_dir,st_meta_dir,meta_key,species,out_f):
    # TODO
    script_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__))) + 'pipelines/'
    r_script_file = f'{script_path}/run_spatalk_lr.R'
    args = [st_dir,st_meta_dir,st_meta_dir,meta_key,species,out_f]
    subprocess.run([rscript_executable, "--vanilla", r_script_file]+ args)