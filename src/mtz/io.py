
class file_reader(object):
  def __init__(self, filename):
    self.filename = filename
    self.stream = open(filename, 'rb')

  def tell(self):
    return self.stream.tell()

  def seek(self, offset, from_what=0):
    self.stream.seek(offset, from_what)

  def close(self):
    self.stream.close()

  def read_constant(self, data):
    filedata = self.stream.read(len(data))
    if not data == filedata:
      raise IOError("Expected constant not encountered; {} != {}".format(filedata, data))

  def read(self, length):
    return self.stream.read(length)
