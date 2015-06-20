#!/usr/bin/env python

import cgi, os
import cgitb; cgitb.enable()

form = cgi.FieldStorage()

# Get filename here.
fileitem = form['filename']

# Test if the file was uploaded
if fileitem.filename :
    fn = os.path.basename(fileitem.filename)
    pic_path = os.path.join('/var/www/ocr/cgi-bin/tmp/', fn)
    out_path = os.path.join('/var/www/ocr/cgi-bin/tmp/', fn)
    open(pic_path, 'wb').write(fileitem.file.read())
    message = 'The file "' + fn + '" was uploaded successfully'

    command = 'tesseract -l eng ' + pic_path + ' ' + out_path
    os.system(command)

    with open(out_path + '.txt', 'r') as f:
        for line in f:
            message += ('<p>' +line + '</p>')

else:
    message = 'No file was uploaded'


print """\
Content-Type: text/html\n
<html>
<body>
   %s
</body>
</html>
""" % (message,)
