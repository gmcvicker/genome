#include <stdio.h>
#include <errno.h>
#include <glib.h>
#include <stdlib.h>

#include "err.h"
#include "util.h"
#include "config.h"
#include "memutil.h"

/* parses a single line of a config file */
static void config_parse_line(Config *conf, char *line) {
  char **array, **tokens;
  char *key, *val;
  int i;

  util_str_lstrip(line);

  /* skip blank or commented lines */
  if(line[0] == '\0' || line[0] == '#') {
    return;
  }

  tokens = g_strsplit(line, "=", 2);

  if(tokens[0] == NULL || tokens[1] == NULL) {
    my_err("%s:%d: expected '=' sign seperating tokens."
	    " line: '%s'", __FILE__, __LINE__, line);
  }

  /* every value is stored as an array and as a string
   * so that it can easily be retrieved as either
   */
  util_str_rstrip(tokens[0]);
  util_str_strip(tokens[1]);

  key = g_strdup(tokens[0]);
  val = g_strdup(tokens[1]);

  g_strfreev(tokens);

  array = g_strsplit(val, ",", 0);
  /* remove leading and trailing whitespace from each array element */
  i = 0;
  while(array[i] != NULL) {
    util_str_strip(array[i]);
    i++;
  }

  /* check that this key is unique */
  if(g_hash_table_lookup(conf->vals, key) != NULL) {
    my_err("%s:%d: configuration key '%s' is defined "
	    "multiple times", __FILE__, __LINE__, key);
  }

  g_hash_table_insert(conf->vals, key, val);
  g_hash_table_insert(conf->val_arrays, key, array);
}


static void free_array(void *ptr) {
  g_strfreev(ptr);
}


/* Reads a configuration file into a hashtable
 * that can then be used to lookup values associated
 * with configuration keys
 */
Config *config_read_file(const char *filename) {
  Config *conf;
  FILE *fh;
  char *line;

  conf = my_new(Config, 1);
  conf->vals = g_hash_table_new_full(g_str_hash, g_str_equal,
				     free, free);

  conf->val_arrays = g_hash_table_new_full(g_str_hash, g_str_equal,
					   NULL, free_array);

  conf->val_long_arrays = g_hash_table_new_full(g_str_hash, g_str_equal,
						NULL, free);

  conf->val_double_arrays = g_hash_table_new_full(g_str_hash, g_str_equal,
						  NULL, free);

  fh = fopen(filename, "r");
  if(fh == NULL) {
    my_err("%s:%d: Could not open file '%s'", __FILE__, __LINE__, filename);
  }
  
  while((line = util_fgets_line(fh)) != NULL) {
    config_parse_line(conf, line);
    my_free(line);
  }
  
  fclose(fh);

  return conf;
}



/**
 * Reads configuration information from command line arguments. The
 * 0th argument is assumed to be the name of the binary, and the 1st
 * argument is path the the configuration file. Remaining arguments
 * are assumed to specify additional config keys formatted as in the
 * same way as the config file: KEY1=VALUE1
 */
Config *config_read_args(const int argc, const char **argv) {
  Config *config;
  char *conf_arg;
  int i;

  if(argc < 2) {
    my_err("%s:%d: expected at least two arguments",
	    __FILE__, __LINE__);
  }
  config = config_read_file(argv[1]);

  for(i = 2; i < argc; i++) {
    conf_arg = g_strdup(argv[i]);
    config_parse_line(config, conf_arg);
    my_free(conf_arg);
  }

  return config;
}


/**
 * Frees the memory allocated for a Config data structure. 
 */
void config_free(Config *conf) {
  g_hash_table_destroy(conf->val_arrays);
  g_hash_table_destroy(conf->val_long_arrays);
  g_hash_table_destroy(conf->val_double_arrays);
  g_hash_table_destroy(conf->vals);
  my_free(conf);

  return;
}



/**
 * returns true if the specified key is defined in the configuration,
 * false otherwise.
 */
int config_has_key(const Config *conf, const char *key) {
  return(g_hash_table_lookup(conf->vals, key) != NULL);
}



/**
 * Gives an error or warning that the config is missing a key
 */
void config_missing_key(const Config *conf, const char *key) {
  my_err("%s:%d: no value associated with config key '%s'",
	 __FILE__, __LINE__, key);
}


/* Returns a the string associated with the
 * provided config key. A warning is issued and NULL returned
 * if the key is not found in the configuration.
 */
char *config_get_str(const Config *conf, const char *key) {
  char *val;

  val = g_hash_table_lookup(conf->vals, key);

  if(val == NULL) {
    config_missing_key(conf, key);
    return NULL;
  }
  
  return val;
}




/* Returns an array of strings associated with the provided
 * configuration key. A warning is issued and NULL returned if the key
 * is not found in the configuration. If the num argument is non-null,
 * the value which it points to is set to be the number of elements in
 * the returned array. The array is also NULL terminated.
 */
char **config_get_str_array(const Config *conf, const char *key, int *num) {
  char **array;

  array = g_hash_table_lookup(conf->val_arrays, key);

  if(array == NULL) {
    config_missing_key(conf, key);

    if(num != NULL) {
      *num = 0;
    }
    return NULL;
  }


  if(num != NULL) {
    *num = 0;
    while(array[*num] != NULL) {
      *num += 1;
    }
  }

  return array;
}





/* Returns a double floating point value associated with the
 * provided configuration key. A warning is issued if the
 * if the key is not found in the configuration or the value
 * does not look like a floating point value.
 */
double config_get_double(const Config *conf, const char *key) {
  char *val;
  double d;

  val = g_hash_table_lookup(conf->vals, key);

  if(val == NULL) {
    config_missing_key(conf, key);
    return 0.0;
  }

  errno = 0;
  d = strtod(val, NULL);
  if(d == 0.0 && errno) {
    my_warn("%s:%d: value associated with key '%s' does not "
	    "look like a double (%s)", __FILE__, __LINE__, key, val);
  }

  return d;
}



/* Returns a long integer value associated with the
 * provided configuration key. A warning is issued if the
 * if the key is not found in the configuration or the value
 * does not look like a long integer value.
 */
long config_get_long(const Config *conf, const char *key) {
  char *val;
  long l;

  val = g_hash_table_lookup(conf->vals, key);

  if(val == NULL) {
    config_missing_key(conf, key);
    return 0;
  }

  errno = 0;
  l = strtol(val, NULL, 10);
  if(errno) {
    my_warn("%s:%d: value associated with key '%s' does not "
	    "look like a long (%s)", __FILE__, __LINE__, key, val);
  }

  return l;
}


/* Returns a long integer value associated with the
 * provided configuration key. A warning is issued if the
 * if the key is not found in the configuration or the value
 * does not look like a long integer value.
 */
int config_get_int(const Config *conf, const char *key) {
  char *val;
  long i;

  val = g_hash_table_lookup(conf->vals, key);

  if(val == NULL) {
    config_missing_key(conf, key);
    return 0;
  }

  errno = 0;
  i = (int)strtol(val, NULL, 10);
  if(errno) {
    my_err("%s:%d: could not parse value '%s' for key '%s' as "
	   "an int (%s)", __FILE__, __LINE__, val, key);
  }

  return i;
}


/* Returns an array of long integers associated with the provided
 * configuration key. A warning is issued and NULL returned if the key
 * is not found in the configuration. If the num argument is non-null,
 * the value which it points to is set to be the number of elements in
 * the returned array.
 */
long *config_get_long_array(Config *conf, const char *key, int *num) {
  char **str_array;
  long *long_array;
  int n;
  int i;

  /* first make sure the array exists as a string array (it is stored
   * as this by default)
   */
  str_array = g_hash_table_lookup(conf->val_arrays, key);
  if(str_array == NULL) {
    config_missing_key(conf, key);

    if(num != NULL) {
      *num = 0;
    }
    return NULL;
  }

  /* count elements in string array */
  n = 0;
  while(str_array[n] != NULL) {
    n++;
  }
  if(num != NULL) {
    *num = n;
  }

  /* check if char array already converted into long data type */
  long_array = g_hash_table_lookup(conf->val_long_arrays, key);
  if(long_array == NULL) {
    /* create the new long array and store it for future use & freeing */
    long_array = my_new(long, n);
    for(i = 0; i < n; i++) {
      errno = 0;
      long_array[i] = strtol(str_array[i], NULL, 10);

      if(errno) {
	my_err("%s:%d: could not parse value '%s' for key '%s' as "
	       "an int (%s)", __FILE__, __LINE__, str_array[i], key);
      }
    }

    g_hash_table_insert(conf->val_long_arrays, g_strdup(key), long_array);
  }
  
  return long_array;
}


/* Returns an array of double floating point numbers associated with
 * the provided configuration key. A warning is issued and NULL
 * returned if the key is not found in the configuration. If the num
 * argument is non-null, the value which it points to is set to be the
 * number of elements in the returned array.
 */
double *config_get_double_array(Config *conf, const char *key, int *num) {
  char **str_array, *next;
  double *dbl_array;
  int n;
  int i;

  /* first make sure the array exists as a string array (it is stored
   * as this by default)
   */
  str_array = g_hash_table_lookup(conf->val_arrays, key);
  if(str_array == NULL) {
    config_missing_key(conf, key);

    if(num != NULL) {
      *num = 0;
    }
    return NULL;
  }

  /* count elements in string array */
  n = 0;
  while(str_array[n] != NULL) {
    n++;
  }
  if(num != NULL) {
    *num = n;
  }

  /* check if char array already converted into long data type */
  dbl_array = g_hash_table_lookup(conf->val_double_arrays, key);
  if(dbl_array == NULL) {
    /* create the new array and store it for future use & freeing */
    dbl_array = my_new(double, n);
    for(i = 0; i < n; i++) {
      dbl_array[i] = strtod(str_array[i], &next);
      if((dbl_array[i] == 0.0) && (next == str_array[i])) {
	my_err("%s:%d: value %s for key %s could not be parsed as a double",
	       __FILE__, __LINE__, str_array[i], key);
      }
    }

    g_hash_table_insert(conf->val_double_arrays, 
			util_str_dup(key), dbl_array);
  }
  
  return dbl_array;
}




/* Returns a boolean value associated with the
 * provided configuration key. A warning is issued if the
 * if the key is not found in the configuration or the value
 * does not look like a boolean value.
 */
int config_get_boolean(const Config *conf, const char *key) {
  char *val;

  val = g_hash_table_lookup(conf->vals, key);

  if(val == NULL) {
    config_missing_key(conf, key);
    return FALSE;
  }

  if((g_strcasecmp(val, "TRUE") == 0) || 
     (g_strcasecmp(val, "T") == 0) ||
     (g_strcasecmp(val, "1") == 0)) {
    return TRUE;
  }
  else if((g_strcasecmp(val, "FALSE") == 0) || 
	  (g_strcasecmp(val, "F") == 0) ||
	  (g_strcasecmp(val, "0") == 0)) {
    return FALSE;
  }

  my_warn("%s:%d: value associated with key '%s' does not "
	  "look like a boolean (%s)", __FILE__, __LINE__, key, val);

  return FALSE;
}
