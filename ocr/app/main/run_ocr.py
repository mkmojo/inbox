import os
import subprocess
from ..models import User, Pic
from .. import db

def run_ocr(user, photo_path, text_path,\
        ocr_engine='tesseract', language = 'eng'):
    """
    photo_path: input/user uploaded file location
    text_path: output/interpreted file location
    !!NOTICE: tesseract appends '.txt' to the end of input filename
    """
    #create directory if not there
    text_dir = os.path.dirname(text_path)
    if not os.path.exists(text_dir):
        os.mkdir(text_dir)

    command = ' '.join([ocr_engine, '-l ' + language, photo_path, text_path])
    try:
        subprocess.call(command.split())
    except :
        return (['Command Fail:\n' + command])

    #save statistics to database
    #need to add in hash later in the future for better security
    pic = Pic(pic_path = photo_path, text_path = text_path + '.txt', owner=user)
    db.session.add(pic)
    db.session.commit()

    message = []
    with open(text_path + '.txt', 'r') as f:
        for line in f:
            line = line.decode('utf-8').strip()
            if line != '':
                message += [line]

    return message
