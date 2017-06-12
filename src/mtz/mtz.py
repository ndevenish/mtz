
import shlex
from .io import file_reader
from collections import namedtuple

HeaderRecord = namedtuple("HeaderRecord", ["keyword", "data"])

def _map_types(string, type_list):
  parts = shlex.split(string)
  assert len(parts) == len(type_list)
  converted = tuple([x(y) for x, y in zip(type_list, parts)])
  return converted


def _parse_record(data):
  keyword = data[:data.find(" ")]
  data = data[len(keyword)+1:].strip()

  if keyword == "VERS":
    return HeaderRecord(keyword, data)
  elif keyword == "TITLE":
    return HeaderRecord(keyword, data)
  elif keyword == "NCOL":
    ints = tuple([int(x) for x in data.split()])
    assert len(ints) == 3
    return HeaderRecord(keyword, ints)
  elif keyword == "CELL":
    print("Waning: Deprecated header entry CELL")
    vals = tuple([float(x) for x in data.split()])
    assert len(vals) == 6
    return HeaderRecord(keyword, vals)
  elif keyword == "SORT":
    vals = tuple([int(x) for x in data.split()])
    assert len(vals) == 5
    return HeaderRecord(keyword, vals)
  elif keyword == "SYMINF":
    parts = shlex.split(data)
    if len(parts) == 7:
      print("Warning: Using undescribed SYMINF format")
      vals = _map_types(data, [int, int, str, int, str, str, str])
    else:  
      vals = _map_types(data, [int, int, str, int, str, str])
    return HeaderRecord(keyword, vals)
  elif keyword == "SYMM":
    return HeaderRecord(keyword, data)
  elif keyword == "RESO":
    vals = tuple([float(x) for x in data.split()])
    assert len(vals) == 2
    return HeaderRecord(keyword, vals)
  elif keyword == "VALM":
    return HeaderRecord(keyword, data)
  elif keyword == "COL" or keyword == "COLUMN":
    vals = _map_types(data, [str, str, float, float, int])
    return HeaderRecord(keyword, vals)
  elif keyword == "NDIF":
    return HeaderRecord(keyword, int(data))
  elif keyword == "PROJECT":
    return HeaderRecord(keyword, _map_types(data, [int, str]))
  elif keyword == "CRYSTAL":
    return HeaderRecord(keyword, _map_types(data, [int, str]))
  elif keyword == "DATASET":
    return HeaderRecord(keyword, _map_types(data, [int, str]))
  elif keyword == "DCELL":
    return HeaderRecord(keyword, _map_types(data, [int, float, float, float, float, float, float]))
  elif keyword == "DWAVEL":
    return HeaderRecord(keyword, _map_types(data, [int, float]))  
  elif keyword == "BATCH":
    return HeaderRecord(keyword, _map_types(data, [int]))  
  else:
    raise IOError("Unrecognised column type " + keyword)


class MTZFile(object):
  def __init__(self, filename):
    self.filename = filename
    self.stream = file_reader(filename)
    self.stream.read_constant(b"MTZ ")
    # Read location of the header record
    header_position = self.stream.read_uint4()
    # Read the number format (after the header record but behaviour is undefined
    # for size other than sizeof(int) == 4)
    number_format = self.stream.read(4)
    assert number_format == b'\x44\x41\x00\x00'

    self.header = self._read_header((header_position-1)*4)

  def _read_header(self, position):
    self.stream.seek(position)
    header_records = []
    # Read a number of 80-character records, ending in END
    while True:
      record = _parse_record(self.stream.read(80).decode("ascii"))
      if record.keyword == "END":
        break
      header_records.append(record)
    return header_records 