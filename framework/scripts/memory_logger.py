#!/usr/bin/env python
import argparse, re, os, sys, subprocess, time, decimal, datetime
from tempfile import TemporaryFile

class Logger():
  def __init__(self):
    self.log = TemporaryFile()
    self.last_position = 0

  def read(self):
    # Setup our cursor position
    cursor_seek = self.log.tell() - 500
    if cursor_seek - 500 <= 0:
      cursor_seek = 0
    # If cursor overlaps previous position, set it to previous position
    if cursor_seek < self.last_position:
      cursor_seek = self.last_position

    # Read log from cursor position
    self.log.seek(cursor_seek)
    output = self.log.read()
    if output != '':
      sys.stdout.write(output)
    self.last_position = self.log.tell()
    return output

  def pstack(self, process_id):
    tmp_proc = subprocess.Popen(['/usr/bin/pstack', str(process_id)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    return tmp_proc.communicate()[0]

class trackPID():
  def __init__(self, process):
    self.args = process
    self.pid = process.track

  def GetMemory(self):
    tmp_proc = subprocess.Popen(['/bin/ps', '-p', str(self.args.track), '-o', 'rss='], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = tmp_proc.communicate()[0]
    if output != '':
      return { 'TOTAL' : int(re.findall(r'(\d+)', str(output))[0]) }
    else:
      return

class sglPID():
  def __init__(self, process):
    self.args = process
    self.logger = Logger()
    self.process = subprocess.Popen(''.join(self.args.run).split(), stdout=self.logger.log, stderr=self.logger.log)
    self.pid = self.process.pid

  def GetMemory(self):
    return_string = {}
    return_string['LOG'] = self.logger.read()
    return_string['PSTACK'] = ''
    if self.args.pstack:
      return_string['PSTACK'] = self.logger.pstack(self.pid)
    tmp_proc = subprocess.Popen(['/bin/ps', '-p', str(self.pid), '-o', 'rss='], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output = tmp_proc.communicate()[0]
    if output != '':
      return_string['TOTAL'] = int(re.findall(r'(\d+)', str(output))[0])
    else:
      return_string['TOTAL'] = 0
    return return_string


class mpiPID():
  def __init__(self, process):
    self.args = process
    self.logger = Logger()
    self.process = subprocess.Popen(''.join(self.args.run).split(), stdout=self.logger.log, stderr=self.logger.log)
    self.pid = self.process.pid

  def _discover_name(self):
    locations = ''.join(self.args.run).split()
    for item in locations:
      if os.path.exists(item):
        return item

  def GetMemory(self):
    one_time = True
    return_string = {}
    pid_dict = {}
    return_string['LOG'] = self.logger.read()
    return_string['PSTACK'] = ''

    command = ['/bin/ps', '-u', str(os.getenv('USER')), '-eo', 'pid,comm']
    tmp_proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    pid_list = tmp_proc.communicate()[0].split('\n')
    for item in pid_list:
      tmp_suffix = item.split('/').pop()
      tmp_item = re.findall(r'(\d+).*' + str(self._discover_name().split('/').pop()), item)
      if len(tmp_item) > 0:
        tmp_proc = subprocess.Popen(['/bin/ps', '-p', str(tmp_item[0]), '-o', 'rss='], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output = tmp_proc.communicate()[0]
        if output != '':
          pid_dict[tmp_item[0]] = int(re.findall(r'(\d+)', str(output))[0])
          if one_time and self.args.pstack:
            return_string['PSTACK'] = self.logger.pstack(str(tmp_item[0]))
            one_time = False
    tmp_mem = 0
    for keyitem in pid_dict:
      tmp_mem += pid_dict[keyitem]
    return_string['TOTAL'] = tmp_mem
    return return_string

class pbsPID():
  def __init__(self, process):
    self.args = process
    self.logfile = Logger()
    self.process = subprocess.Popen(''.join(self.args.run).split(), stdout=self.logger.log, stderr=self.logger.log)
    node_file = open(os.getenv('PBS_NODEFILE'), 'r')
    self.node_list = node_file.read().split()
    node_file.close()
    self.pid = self.process.pid
    self.pid_dict = {}

  def _discover_name(self):
    locations = ''.join(self.args.run).split()
    for item in locations:
      if os.path.exists(item):
        return item

  def GetMemory(self):
    return_string = {}
    return_string['LOG'] = self.logger.read()

    for single_node in self.node_list:
      self.pid_dict[single_node] = {}
      command = ['/usr/bin/ssh', single_node, 'ps', '-u', str(os.getenv('USER')), '-o', 'pid,comm,rss=']
      tmp_proc = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      pid_list = tmp_proc.communicate()[0].split('\n')
      for item in pid_list:
        tmp_suffix = item.split('/').pop()
        tmp_item = re.findall(r'(\d+) .*' + str(self._discover_name().split('/').pop()) +' .* (\d+)', item)
        if len(tmp_item) > 0:
          self.pid_dict[single_node][tmp_item[0][0]] = int(tmp_item[0][1])

    tmp_mem = 0
    for single_node in self.pid_dict:
      for keyitem in self.pid_dict[single_node]:
        tmp_mem += self.pid_dict[single_node][keyitem]
    return_string['TOTAL'] = tmp_mem
    return return_string

class ReadMemoryLog():
  def __init__(self, args):
    history_file = open(args.read, 'r')
    self.memory_list = history_file.read().split('\n')
    self.memory_list.pop()
    history_file.close()
    self.sorted_list = []
    self.mem_list = []
    self.printHistory()

  def printHistory(self):
    RESET  = '\033[0m'
    BOLD   = '\033[1m'
    RED    = '\033[31m'
    GREEN  = '\033[32m'
    CYAN   = '\033[36m'
    YELLOW = '\033[33m'
    last_memory = 0.0
    for timestamp in self.memory_list:
      to = GetTime(eval(timestamp)[0])
      log = eval(timestamp)[2].split('\n')
      pstack = eval(timestamp)[3].split('\n')
      self.mem_list.append(eval(timestamp)[1])
      self.sorted_list.append([str(to.day) + ' ' + str(to.monthname) + ' ' + str(to.hour) + ':' + str(to.minute) + ':' + '{:02.0f}'.format(to.second) + '.' + '{:06.0f}'.format(to.microsecond), eval(timestamp)[1], log, pstack])
    largest_memory = decimal.Decimal(max(self.mem_list))
    percentage_length = decimal.Decimal(self.getTerminalSize()[0]) - decimal.Decimal(len(str(self.sorted_list[0][0]) + ' using: ' + '{:20,.0f}'.format(self.sorted_list[0][1]) + 'K |'))
    print 'Date Stamp' + ' '*int(24) + 'Memory Usage | Percent of MAX memory used: ( ' + str('{:0,.0f}'.format(largest_memory)) + ' K )'
    for item in self.sorted_list:
      if decimal.Decimal(item[1]) == largest_memory:
        percent = '100'
      elif (decimal.Decimal(item[1]) / largest_memory) ==  0:
        percent = '0'
      else:
        percent = str(decimal.Decimal(item[1]) / largest_memory)[2:4] + '.' + str(decimal.Decimal(item[1]) / largest_memory)[4:6]

      tmp_log = ''
      if args.debug:
        for single_log in item[2]:
          if single_log != '':
            tmp_log += ' '*44 + ' | stdout: ' + single_log + '\n'
      if args.pstack:
        for single_pstack in item[3]:
          if single_pstack != '':
            tmp_log += ' '*44 + ' | pstack: ' + single_pstack + '\n'

      if decimal.Decimal(item[1]) == largest_memory:
        tmp_str = item[0] + '{:20,.0f}'.format(item[1]) + ' K |' + '-'*int(percentage_length * (decimal.Decimal(item[1]) / largest_memory)) + RESET + '-| ' + percent + '%\n' + tmp_log
      elif item[1] > last_memory:
        tmp_str = item[0] + '{:20,.0f}'.format(item[1]) + ' K |' + RED + '-'*int(percentage_length * (decimal.Decimal(item[1]) / largest_memory)) + RESET + '| ' + percent + '%\n' + tmp_log
      else:
        tmp_str = item[0] + '{:20,.0f}'.format(item[1]) + ' K |' + GREEN + '-'*int(percentage_length * (decimal.Decimal(item[1]) / largest_memory)) + RESET + '| ' + percent + '%\n' + tmp_log
      last_memory = item[1]
      sys.stdout.write(tmp_str)

    print 'Date Stamp' + ' '*int(22) + 'Memory Usage | Percent of MAX memory used: ( ' + str('{:0,.0f}'.format(largest_memory)) + ' K )'

  def getTerminalSize(self):
    """Quicky to get terminal window size"""
    env = os.environ
    def ioctl_GWINSZ(fd):
      try:
        import fcntl, termios, struct, os
        cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
      except:
        return None
      return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
      try:
        fd = os.open(os.ctermid(), os.O_RDONLY)
        cr = ioctl_GWINSZ(fd)
        os.close(fd)
      except:
        pass
    if not cr:
      try:
        cr = (env['LINES'], env['COLUMNS'])
      except:
        cr = (25, 80)
    return int(cr[1]), int(cr[0])

class GetTime():
  def __init__(self, posix_time=None):
    if posix_time == None:
      self.posix_time = datetime.datetime.now()
    else:
      self.posix_time = datetime.datetime.fromtimestamp(posix_time)
    self.now = float(datetime.datetime.now().strftime('%s.%f'))
    self.microsecond = self.posix_time.microsecond
    self.second = self.posix_time.second
    self.minute = self.posix_time.strftime('%M')
    self.hour = self.posix_time.strftime('%H')
    self.day = self.posix_time.strftime('%d')
    self.month = self.posix_time.strftime('%m')
    self.year = self.posix_time.year
    self.dayname = self.posix_time.strftime('%a')
    self.monthname = self.posix_time.strftime('%b')

class ExportMemoryUsage():
  def __init__(self, args):
    history_file = open(args.export, 'r')
    self.memory_list = history_file.read().split('\n')
    self.memory_list.pop()
    history_file.close()
    self.sorted_list = []
    self.mem_list = []
    self.exportFile()

  def exportFile(self):
    output_file = open(args.export + '.comma_delimited', 'w')
    for timestamp in self.memory_list:
      time_object = GetTime(eval(timestamp)[0])
      self.mem_list.append(eval(timestamp)[1])
      self.sorted_list.append([str(time_object.year) + '-' + str(time_object.month) + '-' + str(time_object.day) + ' ' + str(time_object.hour) + ':' + str(time_object.minute) + ':' + str(time_object.second) + '.' + str(time_object.microsecond), eval(timestamp)[1]])
    for item in self.sorted_list:
      output_file.write(str(item[0]) + ',' + str(item[1]) + '\n')
    output_file.close()
    print 'Comma delimited file saved to: ' + os.getcwd() + '/' + str(args.export) + '.comma_delimited'

def writeMemoryUsage(process):
  file_object = open(process.args.outfile, 'w')
  last_memory = 0
  def _usage(last_memory):
    time_object = GetTime()
    report = process.GetMemory()
    current_usage = [time_object.now, int(report['TOTAL']), report['LOG'], report['PSTACK']]
    if int(current_usage[1]) != int(last_memory):
      last_memory = int(current_usage[1])
      file_object.write(str(current_usage) + '\n')
      file_object.flush()
    # Reset cursor position because the memory has not changed
    else:
      process.logger.last_position = (process.logger.last_position - len(report['LOG']))
    return int(last_memory)
  def _zero():
    log = process.logger.read()
    pstack = ''
    process.logger.log.close()
    time_object = GetTime()
    current_usage = [time_object.now, 0, log, pstack]
    file_object.write(str(current_usage) + '\n')
    file_object.close()
  try:
    while True:
      if not os.kill(int(process.pid), 0):
        if process.args.track is None:
          if process.process.poll() != 0:
            last_memory = _usage(last_memory)
          else:
            _zero()
            time.sleep(1)
            print '\nApplication terminated. Wrote file:', process.args.outfile
            return
        else:
          last_memory = _usage(last_memory)
      else:
        _zero()
        print '\nApplication terminated. Wrote file:', process.args.outfile
        return
      time.sleep(float(process.args.repeat_rate))
  except KeyboardInterrupt:
    _zero()
    print '\nCanceled by user. Wrote file:', process.args.outfile
    return
  except:
    _zero()
    print '\nApplication terminated. Wrote file:', process.args.outfile
    return
#  except:
#    raise

def _verifyARGs(args):
  option_count = 0
  if args.track:
    option_count += 1
  if args.read:
    option_count += 1
  if args.run:
    option_count += 1
  if args.export:
    option_count += 1
  if option_count != 1:
    print 'You must use one of the following: track, read, or run'
    sys.exit(1)
  if args.pbs and args.debug:
    print 'Gathering a stack trace across nodes is not supported.'
    sys.exit(1)
  if args.pstack:
    if not os.path.exists('/usr/bin/pstack'):
      print 'pstack binary not found. You must remove the --debug switch'
      sys.exit(1)
  if args.pbs:
    # There is a limit to how many times we can ssh somewhere (And for good reason)
    # So set the repeat to something that won't flood the network
    if args.repeat_rate <= 5.0:
      args.repeat_rate = 5.0
  if args.run:
    if args.run[0].split(' ')[0] == 'mpiexec' or args.run[0].split(' ')[0] == 'mpirun':
      args.mpi = True
    else:
      args.mpi = False
  return args

def parseARGs(args=None):
  parser = argparse.ArgumentParser(description='Track memory usage')
  parser.add_argument('--track', nargs="?", type=int, help='Track a single specific PID already running\n ')
  parser.add_argument('--read', nargs="?", help='Read a specified memory log file\n ')
  parser.add_argument('--debug', dest='debug', action='store_const', const=True, default=False, help='Display stdout information in output file\n ')
  parser.add_argument('--pstack', dest='pstack', action='store_const', const=True, default=False, help='Save / Display pstack information in output file\n ')
  parser.add_argument('--export', nargs="?", help='Export specified log file to a comma delimited format\n ')
  parser.add_argument('--run', nargs="+", help='Run specified command. You must encapsulate the command in quotes\n ')
  parser.add_argument('--pbs', dest='pbs', action='store_const', const=True, default=False, help='Instruct memory logger to tally all launches on all nodes\n ')
  parser.add_argument('--outfile', nargs="?", default=os.getcwd() + '/usage.log' ,help='Save log to specified file. (Defaults to usage.log)\n ')
  parser.add_argument('--repeat-rate', nargs="?", type=float, default=0.25, help='Indicate the sleep delay in float seconds to check memory usage (default 0.25 seconds)\n ')
  return _verifyARGs(parser.parse_args(args))

if __name__ == '__main__':
  args = parseARGs()
  if args.read:
    ReadMemoryLog(args)
    sys.exit(0)
  if args.run:
    if args.mpi:
      if args.pbs:
        results = pbsPID(args)
      else:
        results = mpiPID(args)
    else:
      results = sglPID(args)
  if args.track:
    results = trackPID(args)
  if args.export:
    ExportMemoryUsage(args)
    sys.exit(0)

  writeMemoryUsage(results)
