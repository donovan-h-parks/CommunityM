#!/usr/bin/env python

###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################


"""
Read configuration file specifying details of all samples.
"""

__author__ = 'Donovan Parks'
__copyright__ = 'Copyright 2013'
__credits__ = ['Donovan Parks']
__license__ = 'GPL3'
__version__ = '1.0.0'
__maintainer__ = 'Donovan Parks'
__email__ = 'donovan.parks@gmail.com'
__status__ = 'Development'

import os, sys, argparse
import ConfigParser

class ReadConfig(object):
  def __init__(self):
    pass

  def readProjectParams(self, configParser):
    projectParams = {}

    for option in configParser.options('project'):
      projectParams[option] = configParser.get('project', option)

    return projectParams

  def readSample(self, configParser, section):
    sampleParams = {}

    for option in configParser.options(section):
      data = configParser.get(section, option)
      if option == 'name':
        sampleParams[option] = data
      else:
        if data != '':
          sampleParams[option] = [x.strip() for x in data.split(',')]
        else:
          sampleParams[option] = []

    return sampleParams

  def validateMothurCompatibility(self, path):
    if '-' in path:
      print '[Error] For compatibility with mothur, paths and files must not contain hypens.'
      print '[Error] The following path and/or file must be renamed: ' + path
      sys.exit()

  def validateConfig(self, configFile, projectParams, allSampleParams, bOutputDirShouldExist):
    # fix output path so it is relative to executing script and ensure it exists
    configPath = configFile[0:configFile.rfind('/')+1]

    projectParams['output_dir'] = configPath + projectParams['output_dir'] + '/'
    projectParams['output_dir'] = projectParams['output_dir'].replace('/./', '/') # make the paths pretty for output
    projectParams['output_dir'] = projectParams['output_dir'].replace('//', '/')
    self.validateMothurCompatibility(projectParams['output_dir'])

    if not bOutputDirShouldExist and os.path.exists(projectParams['output_dir']):
      print '[Error] Output directory already exists.'
      sys.exit()
    elif bOutputDirShouldExist and not os.path.exists(projectParams['output_dir']):
      print '[Error] Output directory does not exists: ' + projectParams['output_dir']
      sys.exit()

    if not bOutputDirShouldExist: # time to create it!
      os.makedirs(projectParams['output_dir'])

    # fix paths so they are relative to the executing script and ensure all files exist
    for sample in allSampleParams:
      pairPaths = []
      for path in allSampleParams[sample]['pairs']:
        pathRelativeToScript = configPath + path
        pathRelativeToScript = pathRelativeToScript.replace('/./', '/') # make the paths pretty for output
        pathRelativeToScript = pathRelativeToScript.replace('//', '/')

        self.validateMothurCompatibility(pathRelativeToScript)

        if not os.path.exists(pathRelativeToScript):
          print '[Error] The following file does not exists: ' + pathRelativeToScript
          sys.exit()

        pairPaths.append(pathRelativeToScript)

      allSampleParams[sample]['pairs'] = pairPaths

      singlePaths = []
      for path in allSampleParams[sample]['singles']:
        pathRelativeToScript = configPath + path
        pathRelativeToScript = pathRelativeToScript.replace('/./', '/') # make the paths pretty for output
        pathRelativeToScript = pathRelativeToScript.replace('//', '/')

        self.validateMothurCompatibility(pathRelativeToScript)

        if not os.path.exists(pathRelativeToScript):
          print '[Error] The following file does not exists: ' + pathRelativeToScript
          sys.exit()

        singlePaths.append(pathRelativeToScript)

      allSampleParams[sample]['singles'] = singlePaths

  def readConfig(self, configFile, outputDirExists):
    configParser = ConfigParser.ConfigParser()
    configParser.readfp(open(configFile))

    projectParams = self.readProjectParams(configParser)

    allSampleParams = {}
    for section in configParser.sections():
      if section[0:6] == 'sample':
        sampleParams = self.readSample(configParser, section)
        allSampleParams[section] = sampleParams

    self.validateConfig(configFile, projectParams, allSampleParams, outputDirExists)

    return projectParams, allSampleParams

  def run(self, configFile):
    projectParams, sampleParams = self.readConfig(configFile)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description="Read configuration file specifying details of all samples.",
                                        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

  parser.add_argument('config_file', help='Project configuration file.')
  args = parser.parse_args()

  readConfig = ReadConfig()
  readConfig.run(args.config_file)
