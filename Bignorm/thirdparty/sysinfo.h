#ifndef PROC_SYSINFO_H
#define PROC_SYSINFO_H
#include <sys/types.h>
#include <dirent.h>
#include "procps.h"

EXTERN_C_BEGIN

extern long smp_num_cpus;          /* number of CPUs */
extern long page_bytes;            /* this architecture's bytes per page */

/* Shmem in 2.6.32+ */
extern unsigned long kb_main_shared;
/* old but still kicking -- the important stuff */
extern unsigned long kb_main_buffers;
extern unsigned long kb_main_cached;
extern unsigned long kb_main_free;
extern unsigned long kb_main_total;
extern unsigned long kb_swap_free;
extern unsigned long kb_swap_total;
/* recently introduced */
extern unsigned long kb_high_free;
extern unsigned long kb_high_total;
extern unsigned long kb_low_free;
extern unsigned long kb_low_total;
extern unsigned long kb_main_available;
/* 2.4.xx era */
extern unsigned long kb_active;
extern unsigned long kb_inact_laundry;  // grrr...
extern unsigned long kb_inact_dirty;
extern unsigned long kb_inact_clean;
extern unsigned long kb_inact_target;
extern unsigned long kb_swap_cached;  /* late 2.4+ */
/* derived values */
extern unsigned long kb_swap_used;
extern unsigned long kb_main_used;
/* 2.5.41+ */
extern unsigned long kb_writeback;
extern unsigned long kb_slab;
extern unsigned long nr_reversemaps;
extern unsigned long kb_committed_as;
extern unsigned long kb_dirty;
extern unsigned long kb_inactive;
extern unsigned long kb_mapped;
extern unsigned long kb_pagetables;

extern void meminfo(void);


extern unsigned long vm_nr_dirty;
extern unsigned long vm_nr_writeback;
extern unsigned long vm_nr_pagecache;
extern unsigned long vm_nr_page_table_pages;
extern unsigned long vm_nr_reverse_maps;
extern unsigned long vm_nr_mapped;
extern unsigned long vm_nr_slab;
extern unsigned long vm_nr_slab_reclaimable;
extern unsigned long vm_nr_slab_unreclaimable;
extern unsigned long vm_nr_active_file;
extern unsigned long vm_nr_inactive_file;
extern unsigned long vm_nr_free_pages;
extern unsigned long vm_pgpgin;
extern unsigned long vm_pgpgout;
extern unsigned long vm_pswpin;
extern unsigned long vm_pswpout;
extern unsigned long vm_pgalloc;
extern unsigned long vm_pgfree;
extern unsigned long vm_pgactivate;
extern unsigned long vm_pgdeactivate;
extern unsigned long vm_pgfault;
extern unsigned long vm_pgmajfault;
extern unsigned long vm_pgscan;
extern unsigned long vm_pgrefill;
extern unsigned long vm_pgsteal;
extern unsigned long vm_kswapd_steal;
extern unsigned long vm_pageoutrun;
extern unsigned long vm_allocstall;

extern void vminfo(void);
extern void cpuinfo (void);

EXTERN_C_END
#endif /* SYSINFO_H */
