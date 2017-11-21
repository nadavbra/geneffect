from __future__ import absolute_import, division, print_function

import os
import traceback

from .util import log

CONFIG_ENV_VARIABLE = 'GENEFFECT_CONFIG_FILE'
DEFAULT_CONFIG_FILE_PATH = '~/.geneffect_config.py'

if CONFIG_ENV_VARIABLE in os.environ:
    config_file_path = os.environ[CONFIG_ENV_VARIABLE]
else:
    config_file_path = DEFAULT_CONFIG_FILE_PATH

config_file_path = os.path.expanduser(config_file_path)

if os.path.exists(config_file_path):
    try:
        raw_config = {}
        exec(open(config_file_path, 'r').read(), raw_config)
    except:
        raise ValueError(('Your configuration file (%s) seems to be invalid. Consider to reset and reconfigure your configuration file by ' + \
                'copying the file default_config.py (which can be found within this module). Caused by error:' + '\n' + '%s') % \
                (config_file_path, traceback.format_exc()))
else:
    raise IOError(('Configuration file doesn\'t exist: %s. You need to create it by copying the file default_config.py (which ' + \
            'can be found within this module).') % config_file_path)
    
class ConfigSetup(object):

    def __init__(self, user_specified_ref_genome):
        self.ref_genome = _get_ref_genome_version(user_specified_ref_genome)
        log('Working with reference genome version %s (user specified: %s).' % (self.ref_genome, user_specified_ref_genome))
        
    def get_variable_value(self, variable_name):
        return _get_variable_value(variable_name, self.ref_genome)
        
    def get_path(self, variable_name):
        return _get_path(variable_name, self.ref_genome)
    
def _get_ref_genome_version(user_specified_ref_genome):
    
    ref_genome_names = _get_raw_variable('REF_GENOME_NAMES')
    
    if user_specified_ref_genome in ref_genome_names:
        return ref_genome_names[user_specified_ref_genome]
    else:
        raise KeyError('Unrecognized reference genome name: %s. Consider modifying REF_GENOME_NAMES in the configuration file.' % user_specified_ref_genome)
    
def _get_path(variable_name, reference_genome_version):
    path = _get_variable_value(variable_name, reference_genome_version)
    assert os.path.exists(path), 'The file path specified in your configuration (%s) doesn\'t exist: %s' % (config_file_path, path)
    return path
    
def _get_variable_value(variable_name, reference_genome_version):
    return _resolve_variable_value(variable_name, _get_raw_variable(variable_name), reference_genome_version)

def _get_raw_variable(variable_name):
    
    global raw_config
    
    try:
        return raw_config[variable_name]
    except KeyError:
        raise KeyError(('Your configuration file (%s) doesn\'t contain a variable %s. Consider to reset and reconfigure your configuration ' + \
                'file by copying the file default_config.py (which can be found within this module).') % (config_file_path, variable_name))
    
def _resolve_variable_value(variable_name, raw_variable, reference_genome_version):
    
    try:
        for option_reference_genome_versions, option_value in raw_variable:
            if reference_genome_version in option_reference_genome_versions:
                return option_value
    except:
        raise ValueError(('The variable %s within your configuration file (%s) seems to be invalid. Caused by error:' + '\n' + '%s') % \
                (variable_name, config_file_path, traceback.format_exc()))
       
    raise KeyError('Reference genome version %s wasn\'t found for the variable %s within your configuration file (%s).' % \
            (reference_genome_version, variable_name, config_file_path))
