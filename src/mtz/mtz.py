

from enum import Enum

import shlex
from .io import file_reader
from collections import namedtuple
import itertools

from pprint import pprint

HeaderRecord = namedtuple("HeaderRecord", ["keyword", "data"])

class DataType(Enum):
  """Data type of represented columns. Maps to:
      H index h,k,l
      J intensity
      F structure amplitude, F
      D anomalous difference
      Q standard deviation of J,F,D or other (but see L and M below)
      G structure amplitude associated with one member of an hkl -h-k-l pair, F(+) or F(-)
      L standard deviation of a column of type G
      K intensity associated with one member of an hkl -h-k-l pair, I(+) or I(-)
      M standard deviation of a column of type K
      E structure amplitude divided by symmetry factor ("epsilon"). Normally scaled as well to give normalised structure factor
      P phase angle in degrees
      W weight (of some sort)
      A phase probability coefficients (Hendrickson/Lattman)
      B BATCH number
      Y M/ISYM, packed partial/reject flag and symmetry number
      I any other integer
      R any other real"""
  Index = "H"
  Intensity = "J"
  StructureAmplitude = "F"
  AnomalousDifference = "D"
  StandardDeviation  = "Q"
  MemberStructureAmplitude = "G"
  MSAStandardDeviation = "L"
  IndexMemberIntensity = "K"
  IMIStandardDeviation = "M"
  Epsilon = "E"
  PhaseAngle = "P"
  Weight = "W"
  PhaseProbabilityCoefficient = "A"
  BatchNumber = "B"
  MISYM = "Y"
  Integer = "I"
  Real = "R"

def only(iterator):
  generator = (x for x in iterator)
  try:
    value = next(generator)
  except StopIteration:
    raise ValueError("No items in generator to extract")
  try:
    next(generator)
  except StopIteration:
    pass
  else:
    raise ValueError("More than one item in generator")
  return value

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

class MTZFileError(IOError):
  pass

class InconsistentHeaderError(MTZFileError):
  pass

def _parse_batch(stream):
  """Parse a batch, or fail if we've reached the end of them"""
  header = _read_record(stream, expected="BH")

  serial, nwords, nintegers, nreals = header.data
  title = _read_record(stream, expected="TITLE")

  # Read the first three ints and match them
  assert nwords     == stream.read_uint4()
  assert nintegers  == stream.read_uint4()
  assert nreals     == stream.read_uint4()
  assert nreals+nintegers == nwords
  integers = [stream.read_uint4() for x in range(nintegers-3)]
  reals    = [stream.read_float4() for x in range(nreals)]
  
  BHCH = _read_record(stream, expected="BHCH")
  
  return Batch(serial, title.data, (integers, reals), BHCH.data)

def _read_record(stream, expected=None):
  record = _parse_record(stream.read(80).decode("ascii"))
  if expected is not None and record.keyword != expected:
    raise MTZFileError("Unexpected result looking for {} record".format(expected))
  return record
  
def _parse_header(stream):
  header_records = []
  # Read a number of 80-character records, ending in END
  while True:
    record = _read_record(stream)
    if record.keyword == "END":
      break
    header_records.append(record)


  ncols, nrefl, nbatch = only([x for x in header_records if x.keyword == "NCOL"]).data

  # Read the next entry
  history = []
  batches = []
  while True:
    postheader = _read_record(stream)
    
    if postheader.keyword == "MTZHIST":
      history = [stream.read(80).decode("ascii").strip() for x in range(postheader.data)]
    elif postheader.keyword == "MTZBATS":
      batches = [_parse_batch(stream) for i in range(nbatch)]
    if postheader.keyword == "MTZENDOFHEADERS":
      break

  # Consolidate header records with batch records and ensure they match
  header_batch_entries = list(itertools.chain(*[x.data for x in header_records if x.keyword == "BATCH"]))
  if len(header_batch_entries) != len(batches):
    raise InconsistentHeaderError("Number of batches does not match declared header batch index count")
  # Validate the header BATCH serial numbers match the batches BH serial
  for entry, batch in zip(header_batch_entries, batches):
    if not entry == batch.serial:
      raise InconsistentHeaderError("Batch serial order does not match header batch serial order")
  # At this point, we can discard the header batch entries
  header_records = [x for x in header_records if x.keyword != "BATCH"]

  if ncols != len([x for x in header_records if x.keyword in ("COL", "COLUMN")]):
    raise InconsistentHeaderError("Column entries in header do not match column count")

  return Header(header_records, history=history, batches=batches) 


def _parse_record(raw_data):
  keyword = raw_data[:raw_data.find(" ")]
  data = raw_data[len(keyword)+1:]

  if keyword == "VERS":
    return HeaderRecord(keyword, data.strip())
  elif keyword == "TITLE":
    return HeaderRecord(keyword, data.strip())
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
      vals = _map_types(data, [int, int, unicode, int, unicode, unicode, unicode])
    else:  
      vals = _map_types(data, [int, int, unicode, int, unicode, unicode])
    return HeaderRecord(keyword, vals)
  elif keyword == "SYMM":
    return HeaderRecord(keyword, data.strip())
  elif keyword == "RESO":
    vals = tuple([float(x) for x in data.split()])
    assert len(vals) == 2
    return HeaderRecord(keyword, vals)
  elif keyword == "VALM":
    return HeaderRecord(keyword, data.strip())
  elif keyword == "COL" or keyword == "COLUMN":
    vals = _map_types(data, [unicode, unicode, float, float, int])
    return HeaderRecord(keyword, vals)
  elif keyword == "NDIF":
    return HeaderRecord(keyword, int(data))
  elif keyword == "PROJECT":
    return HeaderRecord(keyword, _map_types(data, [int, unicode]))
  elif keyword == "CRYSTAL":
    return HeaderRecord(keyword, _map_types(data, [int, unicode]))
  elif keyword == "DATASET":
    return HeaderRecord(keyword, _map_types(data, [int, unicode]))
  elif keyword == "DCELL":
    return HeaderRecord(keyword, _map_types(data, [int, float, float, float, float, float, float]))
  elif keyword == "DWAVEL":
    return HeaderRecord(keyword, _map_types(data, [int, float]))  
  elif keyword == "BATCH":
    return HeaderRecord(keyword, tuple([int(x) for x in data.split()]))  
  elif keyword == "COLSRC":
    vals = data[:30].strip(), data[31:67].strip(), int(data[66:])
    return HeaderRecord(keyword, vals)
  elif keyword == "COLGRP":
    # COLGRP %-30s %-30s %-4s %1X %4d"
    vals = split_length(data, [31, 31, 5, 2, 4])
    hex_int = lambda x: int(x, 16)
    return HeaderRecord(keyword, _map_types(vals, [unicode, unicode, unicode, hex_int, int]))
  elif keyword == "END":
    return HeaderRecord(keyword, None)
  elif keyword == "MTZHIST":
    return HeaderRecord(keyword, int(data.strip()))
  elif keyword == "MTZENDOFHEADERS":
    return HeaderRecord(keyword, None)
  elif keyword == "MTZBATS":
    return HeaderRecord(keyword, None)
  elif keyword == "BH":
    return HeaderRecord(keyword, _map_types(data, [int, int, int, int]))
  elif keyword == "BHCH":
    vals = [x for x in split_length(data, [9,9,9]) if x]
    return HeaderRecord(keyword, vals)
  else:
    raise IOError("Unrecognised column type " + keyword)


class Batch(object):
  def __init__(self, serial, title, data, bhch):
    self.serial = serial
    self.title = title
    self.data = data
    self.bhch = bhch

class Header(object):
  def __init__(self, raw_records, history=[], batches=[]):
    self.raw_records = raw_records
    self.history = history
    self.batches = batches



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

    print ("\n" + filename)
    pprint(self.header.raw_records)