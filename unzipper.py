import zipfile as zip
import os
from os.path import join, splitext
from glob import glob, iglob


def unzip_directory(zip_dir, target_dir):
    CWD = os.getcwd()
    zip_path = CWD + '/' + zip_dir
    target_path = CWD + '/' + target_dir
    if target_dir != zip_dir:
        os.mkdir(target_dir)

    for zf in iglob(join(zip_path, '*.zip')):

        zip_handler = zip.ZipFile(zf, 'r')
        zip_handler.extractall(target_path)
