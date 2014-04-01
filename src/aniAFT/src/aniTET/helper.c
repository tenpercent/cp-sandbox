#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>
#include "helper.h"

void *libaft_malloc(long size) {
    void *p = malloc(size);
    if (!p) {
	perror("libaft: malloc");
	return 0;
    }
    return p;
}

void *libaft_realloc(void *p, long size) {
    p = realloc(p, size);
    if (!p) {
	perror("libaft: realloc");
	return 0;
    }
    return p;
}

void libaft_free(void*p) {
    free(p);
}

void libaft_2d_warn(char *str, ...) {
    va_list ArgList;
    va_start(ArgList, str);
    fprintf(stderr, "\nlibaft(2d) warning: ");
    vfprintf(stderr, str, ArgList);
    fprintf(stderr, "\n");
}

void libaft_2d_stop(char *str, ...) {
    va_list ArgList;
    va_start(ArgList, str);
    fprintf(stderr, "\nlibaft(2d) error\n~~~~~~~~~~~~~~~~\n");
    vfprintf(stderr, str, ArgList);
    fprintf(stderr, "\n");
//    exit(1);
}

void libaft_3d_warn(char *str, ...) {
    va_list ArgList;
    va_start(ArgList, str);
    fprintf(stderr, "\nlibaft(3d) warning: ");
    vfprintf(stderr, str, ArgList);
    fprintf(stderr, "\n");
}

void libaft_3d_stop(char *str, ...) {
    va_list ArgList;
    va_start(ArgList, str);
    fprintf(stderr, "\nlibaft(3d) error\n~~~~~~~~~~~~~~~~\n");
    vfprintf(stderr, str, ArgList);
    fprintf(stderr, "\n");
//    exit(1);
}

