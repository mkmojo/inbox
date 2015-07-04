from flask.ext.wtf import Form
from wtforms import StringField, SubmitField, FileField
from wtforms.validators import Required

class UploadForm(Form):
    photo = FileField('Which receipt would you love to upload ?', validators=[Required()])
    submit = SubmitField('Upload')
