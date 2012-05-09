#ifndef __CONFIG_H__
#define __CONFIG_H__


#include <glib.h>


typedef struct {
  int missing_key_action;
  GHashTable *vals;
  GHashTable *val_arrays;
  GHashTable *val_long_arrays;
  GHashTable *val_double_arrays;
} Config;


Config *config_read_file(const char *filename);
Config *config_read_args(const int argc, const char **argv);
int config_has_key(const Config *conf, const char *key);
char  *config_get_str(const Config *conf, const char *key);
char **config_get_str_array(const Config *conf, const char *key, int *sz);
double config_get_double(const Config *conf, const char *key);
double *config_get_double_array(Config *conf, const char *key, int *sz);
long   config_get_long(const Config *conf, const char *key);
int   config_get_int(const Config *conf, const char *key);
long  *config_get_long_array(Config *conf, const char *key, int *sz);
int   config_get_boolean(const Config *conf, const char *key);
void config_free(Config *conf);

#endif
