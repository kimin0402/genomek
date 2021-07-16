import subprocess

def run_blat_server(port:int=2882, 
                    blat_path:str='/home/users/kimin/tools/blat/gfServer',
                    reference_path:str="/home/users/kimin/projects/00_Reference/hg19.2bit"):

    cmd = f"{blat_path} start 127.0.0.1 {port} -stepSize=5 -log=/home/users/kimin/projects/00_Reference/blat.log {reference_path}"

    subprocess.run(cmd, shell=True)

    return None