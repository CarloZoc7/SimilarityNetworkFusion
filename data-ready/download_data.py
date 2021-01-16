import os
import requests
import zipfile
import shutil

if __name__ == '__main__':
    to_file = 'https://www.dropbox.com/s/c8uixr6ysj1kj3d/data-ready.zip?dl=0'
    url = "https://uc18aedf8001d4e4ebada39bbfc4.dl.dropboxusercontent.com/cd/0/get/BHEFnmRw_MJ_OQysqzG_vdtOlM6_I0GWFoG48DJgfRmSU8NTreDl7JA2vEn7PsRdOG61NzzHEP90Do00m8ClI5mC-MghaUT-V1UWXfn5pvaEtZs8SE2JPoOLX8yBh-Eldso/file?dl=1#"
    name_file = './data-ready.zip'

    print('Downloading the file...', end='')
    
    with requests.get(url, stream=True) as r:
        with open(name_file, 'wb') as f:
            shutil.copyfileobj(r.raw, f)
    print('DONE')
    print('Extracting the file...', end='')
    with zipfile.ZipFile(name_file, 'r') as zip_ref:
        zip_ref.extractall('./')
    os.remove(name_file)
    print('END')
    print('Bye!')