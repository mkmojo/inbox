import os
import subprocess
from ..models import User, Pic
from .. import db
from manage import app
from time import gmtime, strftime
from flask import url_for
from werkzeug import secure_filename


def run_ocr(user, form, \
        ocr_engine='tesseract', language = 'eng'):
    """
    !!NOTICE: tesseract appends '.txt' to the end of input filename
    """

    #prefix timestamp to filename
    filename = strftime("%Y-%m-%d-%H:%M:%S_", gmtime()) + \
            secure_filename(form.photo.data.filename)

    #relative path to static folder
    relative_pic_path = os.path.join('uploads/', str(user.id),  filename)
    pic_path = os.path.join(app.config['STATIC_FOLDER'], relative_pic_path)

    relative_text_path = os.path.join('output/', str(user.id),  filename)
    text_path = os.path.join(app.config['STATIC_FOLDER'], relative_text_path)

    pic_dir_name = os.path.dirname(pic_path)
    try:
        if not os.path.exists(pic_dir_name):
            os.mkdir(pic_dir_name)
        form.photo.data.save(pic_path)
    except:
        print ("failed to save picture to: " + pic_path)

    #create directory if not there
    text_dir = os.path.dirname(text_path)
    if not os.path.exists(text_dir):
        os.mkdir(text_dir)

    command = ' '.join([ocr_engine, '-l ' + language, pic_path, text_path])
    try:
        subprocess.call(command.split())
    except :
        return (['Command Fail:\n' + command])

    #save statistics to database
    #need to add in hash later in the future for better security
    pic = Pic(pic_path = relative_pic_path, text_path = relative_text_path + \
            '.txt', owner=user)
    db.session.add(pic)
    db.session.commit()

    message = []
    with open(text_path + '.txt', 'r') as f:
        for line in f:
            line = line.decode('utf-8').strip()
            if line != '':
                message += [line]

    return message
