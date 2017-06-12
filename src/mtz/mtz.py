
import shlex
from .io import file_reader
from collections import namedtuple

HeaderRecord = namedtuple("HeaderRecord", ["keyword", "data"])

def _map_types(string, type_list):
  if isinstance(string, basestring):
    parts = shlex.split(string)
  else:
    parts = string
  assert len(parts) == len(type_list)
  converted = tuple([x(y) for x, y in zip(type_list, parts)])
  return converted

def split_length(string, lengths):
  """Split a string based on lengths.
  Designed to read c-scanf-based output"""
  parts = []
  pos = 0
  for length in lengths:
    parts.append(string[pos:pos+length].strip())
    pos += length
  return parts

def _parse_header(stream):
  header_records = []
  # Read a number of 80-character records, ending in END
  while True:
    record = _parse_record(stream.read(80).decode("ascii"))
    if record.keyword == "END":
      break
    header_records.append(record)

  # Read the next entry
  postheader = _parse_record(stream.read(80).decode("ascii"))
  if postheader.keyword == "MTZHIST":
    history = [stream.read(80).decode("ascii").strip() for x in range(postheader.data)]

  posthistory = _parse_record(stream.read(80).decode("ascii"))
  assert posthistory.keyword == "MTZENDOFHEADERS"

  return MTZHeader(header_records, history=history) 


def _parse_record(raw_data):
  keyword = raw_data[:raw_data.find(" ")]
  data = raw_data[len(keyword)+1:]

  if keyword == "VERS":
    return HeaderRecord(keyword, data)
  elif keyword == "TITLE":
    return HeaderRecord(keyword, data)
  elif keyword == "NCOL":
    ints = tuple([int(x) for x in data.split()])
    assert len(ints) == 3
    return HeaderRecord(keyword, ints)
  elif keyword == "CELL":
    # print("Warning: Deprecated header entry CELL")
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
 #       sprintf(hdrrec,"SYMINF %3d %2d %c %5d %22s %5s %c",mtz->mtzsymm.nsym,
 # 2713           mtz->mtzsymm.nsymp,mtz->mtzsymm.symtyp,mtz->mtzsymm.spcgrp,spgname,
 # 2714           mtz->mtzsymm.pgname,mtz->mtzsymm.spg_confidence);
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
  elif keyword == "COLSRC":
    vals = data[:30].strip(), data[31:67].strip(), int(data[66:])
    return HeaderRecord(keyword, vals)
  elif keyword == "COLGRP":
    # COLGRP %-30s %-30s %-4s %1X %4d"
    vals = split_length(data, [31, 31, 5, 2, 4])
    hex_int = lambda x: int(x, 16)
    return HeaderRecord(keyword, _map_types(vals, [str, str, str, hex_int, int]))
  elif keyword == "END":
    return HeaderRecord(keyword, None)
  elif keyword == "MTZHIST":
    return HeaderRecord(keyword, int(data.strip()))
  elif keyword == "MTZENDOFHEADERS":
    return HeaderRecord(keyword, None)
  else:
    raise IOError("Unrecognised column type " + keyword)


class MTZHeader(object):
  def __init__(self, raw_records, history=[]):
    self.raw_records = raw_records
    self.history = history

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

    self.stream.seek((header_position-1)*4)
    self.header = _parse_header(self.stream)