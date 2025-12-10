/*
 * File for parsing top-level /proc entities.
 * Copyright (C) 1992-1998 by Michael K. Johnson, johnsonm@redhat.com
 * Copyright 1998-2003 Albert Cahalan
 * June 2003, Fabian Frederick, disk and slab info
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */


// unneeded parts deleted

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <locale.h>

#include <unistd.h>
#include <fcntl.h>
#include "alloc.h"
// #include "version.h"
#include "sysinfo.h" /* include self to verify prototypes */

#ifndef HZ
#include <netinet/in.h>  /* htons */
#endif

long smp_num_cpus;     /* number of CPUs */
long page_bytes;       /* this architecture's page size */

#define BAD_OPEN_MESSAGE					\
"Error: /proc must be mounted\n"				\
"  To mount /proc at boot you need an /etc/fstab line like:\n"	\
"      proc   /proc   proc    defaults\n"			\
"  In the meantime, run \"mount proc /proc -t proc\"\n"

#define STAT_FILE    "/proc/stat"
static int stat_fd = -1;
#define UPTIME_FILE  "/proc/uptime"
static int uptime_fd = -1;
#define LOADAVG_FILE "/proc/loadavg"
static int loadavg_fd = -1;
#define MEMINFO_FILE "/proc/meminfo"
static int meminfo_fd = -1;
#define VMINFO_FILE "/proc/vmstat"
static int vminfo_fd = -1;
#define VM_MIN_FREE_FILE "/proc/sys/vm/min_free_kbytes"
static int vm_min_free_fd = -1;

// As of 2.6.24 /proc/meminfo seems to need 888 on 64-bit,
// and would need 1258 if the obsolete fields were there.
// As of 3.13 /proc/vmstat needs 2623,
// and /proc/stat needs 3076.
static char buf[8192];

/* This macro opens filename only if necessary and seeks to 0 so
 * that successive calls to the functions are more efficient.
 * It also reads the current contents of the file into the global buf.
 */
#define FILE_TO_BUF(filename, fd) do{				\
    static int local_n;						\
    if (fd == -1 && (fd = open(filename, O_RDONLY)) == -1) {	\
	fputs(BAD_OPEN_MESSAGE, stderr);			\
	fflush(NULL);						\
	_exit(102);						\
    }								\
    lseek(fd, 0L, SEEK_SET);					\
    if ((local_n = read(fd, buf, sizeof buf - 1)) < 0) {	\
	perror(filename);					\
	fflush(NULL);						\
	_exit(103);						\
    }								\
    buf[local_n] = '\0';					\
}while(0)

/* evals 'x' twice */
#define SET_IF_DESIRED(x,y) do{  if(x) *(x) = (y); }while(0)

/* return minimum of two values */
#define MIN(x,y) ((x) < (y) ? (x) : (y))


/***********************************************************************/
/*
 * Copyright 1999 by Albert Cahalan; all rights reserved.
 * This file may be used subject to the terms and conditions of the
 * GNU Library General Public License Version 2, or any later version
 * at your option, as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Library General Public License for more details.
 */

typedef struct mem_table_struct {
  const char *name;     /* memory type name */
  unsigned long *slot; /* slot in return struct */
} mem_table_struct;

static int compare_mem_table_structs(const void *a, const void *b){
  return strcmp(((const mem_table_struct*)a)->name,((const mem_table_struct*)b)->name);
}


/* Shmem in 2.6.32+ */
unsigned long kb_main_shared;
/* old but still kicking -- the important stuff */
static unsigned long kb_page_cache;
unsigned long kb_main_buffers;
unsigned long kb_main_free;
unsigned long kb_main_total;
unsigned long kb_swap_free;
unsigned long kb_swap_total;
/* recently introduced */
unsigned long kb_high_free;
unsigned long kb_high_total;
unsigned long kb_low_free;
unsigned long kb_low_total;
unsigned long kb_main_available;
/* 2.4.xx era */
unsigned long kb_active;
unsigned long kb_inact_laundry;
unsigned long kb_inact_dirty;
unsigned long kb_inact_clean;
unsigned long kb_inact_target;
unsigned long kb_swap_cached;  /* late 2.4 and 2.6+ only */
/* derived values */
unsigned long kb_main_cached;
unsigned long kb_swap_used;
unsigned long kb_main_used;
/* 2.5.41+ */
unsigned long kb_writeback;
unsigned long kb_slab;
unsigned long nr_reversemaps;
unsigned long kb_committed_as;
unsigned long kb_dirty;
unsigned long kb_inactive;
unsigned long kb_mapped;
unsigned long kb_pagetables;
// seen on a 2.6.x kernel:
static unsigned long kb_vmalloc_chunk;
static unsigned long kb_vmalloc_total;
static unsigned long kb_vmalloc_used;
// seen on 2.6.24-rc6-git12
static unsigned long kb_anon_pages;
static unsigned long kb_bounce;
static unsigned long kb_commit_limit;
static unsigned long kb_nfs_unstable;
// seen on 2.6.18
static unsigned long kb_min_free;
// 2.6.19+
static unsigned long kb_slab_reclaimable;
static unsigned long kb_slab_unreclaimable;
// 2.6.27+
static unsigned long kb_active_file;
static unsigned long kb_inactive_file;


void meminfo(void){
  char namebuf[32]; /* big enough to hold any row name */
  mem_table_struct findme = { namebuf, NULL};
  mem_table_struct *found;
  char *head;
  char *tail;
  static const mem_table_struct mem_table[] = {
  {"Active",       &kb_active},       // important
  {"Active(file)", &kb_active_file},
  {"AnonPages",    &kb_anon_pages},
  {"Bounce",       &kb_bounce},
  {"Buffers",      &kb_main_buffers}, // important
  {"Cached",       &kb_page_cache},  // important
  {"CommitLimit",  &kb_commit_limit},
  {"Committed_AS", &kb_committed_as},
  {"Dirty",        &kb_dirty},        // kB version of vmstat nr_dirty
  {"HighFree",     &kb_high_free},
  {"HighTotal",    &kb_high_total},
  {"Inact_clean",  &kb_inact_clean},
  {"Inact_dirty",  &kb_inact_dirty},
  {"Inact_laundry",&kb_inact_laundry},
  {"Inact_target", &kb_inact_target},
  {"Inactive",     &kb_inactive},     // important
  {"Inactive(file)",&kb_inactive_file},
  {"LowFree",      &kb_low_free},
  {"LowTotal",     &kb_low_total},
  {"Mapped",       &kb_mapped},       // kB version of vmstat nr_mapped
  {"MemAvailable", &kb_main_available}, // important
  {"MemFree",      &kb_main_free},    // important
  {"MemTotal",     &kb_main_total},   // important
  {"NFS_Unstable", &kb_nfs_unstable},
  {"PageTables",   &kb_pagetables},   // kB version of vmstat nr_page_table_pages
  {"ReverseMaps",  &nr_reversemaps},  // same as vmstat nr_page_table_pages
  {"SReclaimable", &kb_slab_reclaimable}, // "slab reclaimable" (dentry and inode structures)
  {"SUnreclaim",   &kb_slab_unreclaimable},
  {"Shmem",        &kb_main_shared},  // kernel 2.6.32 and later
  {"Slab",         &kb_slab},         // kB version of vmstat nr_slab
  {"SwapCached",   &kb_swap_cached},
  {"SwapFree",     &kb_swap_free},    // important
  {"SwapTotal",    &kb_swap_total},   // important
  {"VmallocChunk", &kb_vmalloc_chunk},
  {"VmallocTotal", &kb_vmalloc_total},
  {"VmallocUsed",  &kb_vmalloc_used},
  {"Writeback",    &kb_writeback},    // kB version of vmstat nr_writeback
  };
  const int mem_table_count = sizeof(mem_table)/sizeof(mem_table_struct);
  unsigned long watermark_low;
  signed long mem_available;

  FILE_TO_BUF(MEMINFO_FILE,meminfo_fd);

  kb_inactive = ~0UL;
  kb_low_total = kb_main_available = 0;

  head = buf;
  for(;;){
    tail = strchr(head, ':');
    if(!tail) break;
    *tail = '\0';
    if(strlen(head) >= sizeof(namebuf)){
      head = tail+1;
      goto nextline;
    }
    strcpy(namebuf,head);
    found = bsearch(&findme, mem_table, mem_table_count,
        sizeof(mem_table_struct), compare_mem_table_structs
    );
    head = tail+1;
    if(!found) goto nextline;
    *(found->slot) = (unsigned long)strtoull(head,&tail,10);
nextline:
    tail = strchr(head, '\n');
    if(!tail) break;
    head = tail+1;
  }
  if(!kb_low_total){  /* low==main except with large-memory support */
    kb_low_total = kb_main_total;
    kb_low_free  = kb_main_free;
  }
  if(kb_inactive==~0UL){
    kb_inactive = kb_inact_dirty + kb_inact_clean + kb_inact_laundry;
  }
  kb_main_cached = kb_page_cache + kb_slab;
  kb_swap_used = kb_swap_total - kb_swap_free;
  kb_main_used = kb_main_total - kb_main_free - kb_main_cached - kb_main_buffers;

  /* zero? might need fallback for 2.6.27 <= kernel <? 3.14 */
  if (!kb_main_available) {
      FILE_TO_BUF(VM_MIN_FREE_FILE, vm_min_free_fd);
      kb_min_free = (unsigned long) strtoull(buf,&tail,10);

      watermark_low = kb_min_free * 5 / 4; /* should be equal to sum of all 'low' fields in /proc/zoneinfo */

      mem_available = (signed long)kb_main_free - watermark_low
      + kb_inactive_file + kb_active_file - MIN((kb_inactive_file + kb_active_file) / 2, watermark_low)
      + kb_slab_reclaimable - MIN(kb_slab_reclaimable / 2, watermark_low);

      if (mem_available < 0) mem_available = 0;
      kb_main_available = (unsigned long)mem_available;
  }
}

/*****************************************************************/

/* read /proc/vminfo only for 2.5.41 and above */

typedef struct vm_table_struct {
  const char *name;     /* VM statistic name */
  unsigned long *slot;       /* slot in return struct */
} vm_table_struct;

static int compare_vm_table_structs(const void *a, const void *b){
  return strcmp(((const vm_table_struct*)a)->name,((const vm_table_struct*)b)->name);
}

// see include/linux/page-flags.h and mm/page_alloc.c
unsigned long vm_nr_dirty;           // dirty writable pages
unsigned long vm_nr_writeback;       // pages under writeback
unsigned long vm_nr_pagecache;       // pages in pagecache -- gone in 2.5.66+ kernels
unsigned long vm_nr_page_table_pages;// pages used for pagetables
unsigned long vm_nr_reverse_maps;    // includes PageDirect
unsigned long vm_nr_mapped;          // mapped into pagetables
unsigned long vm_nr_slab;            // in slab
unsigned long vm_nr_slab_reclaimable;  // 2.6.19+ kernels
unsigned long vm_nr_slab_unreclaimable;// 2.6.19+ kernels
unsigned long vm_nr_active_file;       // 2.6.27+ kernels
unsigned long vm_nr_inactive_file;     // 2.6.27+ kernels
unsigned long vm_nr_free_pages;        // 2.6.21+ kernels
unsigned long vm_pgpgin;             // kB disk reads  (same as 1st num on /proc/stat page line)
unsigned long vm_pgpgout;            // kB disk writes (same as 2nd num on /proc/stat page line)
unsigned long vm_pswpin;             // swap reads     (same as 1st num on /proc/stat swap line)
unsigned long vm_pswpout;            // swap writes    (same as 2nd num on /proc/stat swap line)
unsigned long vm_pgalloc;            // page allocations
unsigned long vm_pgfree;             // page freeings
unsigned long vm_pgactivate;         // pages moved inactive -> active
unsigned long vm_pgdeactivate;       // pages moved active -> inactive
unsigned long vm_pgfault;           // total faults (major+minor)
unsigned long vm_pgmajfault;       // major faults
unsigned long vm_pgscan;          // pages scanned by page reclaim
unsigned long vm_pgrefill;       // inspected by refill_inactive_zone
unsigned long vm_pgsteal;       // total pages reclaimed
unsigned long vm_kswapd_steal; // pages reclaimed by kswapd
// next 3 as defined by the 2.5.52 kernel
unsigned long vm_pageoutrun;  // times kswapd ran page reclaim
unsigned long vm_allocstall; // times a page allocator ran direct reclaim
unsigned long vm_pgrotated; // pages rotated to the tail of the LRU for immediate reclaim
// seen on a 2.6.8-rc1 kernel, apparently replacing old fields
static unsigned long vm_pgalloc_dma;          //
static unsigned long vm_pgalloc_high;         //
static unsigned long vm_pgalloc_normal;       //
static unsigned long vm_pgrefill_dma;         //
static unsigned long vm_pgrefill_high;        //
static unsigned long vm_pgrefill_normal;      //
static unsigned long vm_pgscan_direct_dma;    //
static unsigned long vm_pgscan_direct_high;   //
static unsigned long vm_pgscan_direct_normal; //
static unsigned long vm_pgscan_kswapd_dma;    //
static unsigned long vm_pgscan_kswapd_high;   //
static unsigned long vm_pgscan_kswapd_normal; //
static unsigned long vm_pgsteal_dma;          //
static unsigned long vm_pgsteal_high;         //
static unsigned long vm_pgsteal_normal;       //
// seen on a 2.6.8-rc1 kernel
static unsigned long vm_kswapd_inodesteal;    //
static unsigned long vm_nr_unstable;          //
static unsigned long vm_pginodesteal;         //
static unsigned long vm_slabs_scanned;        //

void vminfo(void){
  char namebuf[32]; /* big enough to hold any row name */
  vm_table_struct findme = { namebuf, NULL};
  vm_table_struct *found;
  char *head;
  char *tail;
  static const vm_table_struct vm_table[] = {
  {"allocstall",          &vm_allocstall},
  {"kswapd_inodesteal",   &vm_kswapd_inodesteal},
  {"kswapd_steal",        &vm_kswapd_steal},
  {"nr_active_file",      &vm_nr_active_file},     // 2.6.27+ kernels
  {"nr_dirty",            &vm_nr_dirty},           // page version of meminfo Dirty
  {"nr_free_pages",       &vm_nr_free_pages},      // 2.6.21+ kernels
  {"nr_inactive_file",    &vm_nr_inactive_file},   // 2.6.27+ kernels
  {"nr_mapped",           &vm_nr_mapped},          // page version of meminfo Mapped
  {"nr_page_table_pages", &vm_nr_page_table_pages},// same as meminfo PageTables
  {"nr_pagecache",        &vm_nr_pagecache},       // gone in 2.5.66+ kernels
  {"nr_reverse_maps",     &vm_nr_reverse_maps},    // page version of meminfo ReverseMaps GONE
  {"nr_slab",             &vm_nr_slab},            // page version of meminfo Slab (gone in 2.6.19+)
  {"nr_slab_reclaimable", &vm_nr_slab_reclaimable},// 2.6.19+ kernels
 {"nr_slab_unreclaimable",&vm_nr_slab_unreclaimable},// 2.6.19+ kernels
  {"nr_unstable",         &vm_nr_unstable},
  {"nr_writeback",        &vm_nr_writeback},       // page version of meminfo Writeback
  {"pageoutrun",          &vm_pageoutrun},
  {"pgactivate",          &vm_pgactivate},
  {"pgalloc",             &vm_pgalloc},  // GONE (now separate dma,high,normal)
  {"pgalloc_dma",         &vm_pgalloc_dma},
  {"pgalloc_high",        &vm_pgalloc_high},
  {"pgalloc_normal",      &vm_pgalloc_normal},
  {"pgdeactivate",        &vm_pgdeactivate},
  {"pgfault",             &vm_pgfault},
  {"pgfree",              &vm_pgfree},
  {"pginodesteal",        &vm_pginodesteal},
  {"pgmajfault",          &vm_pgmajfault},
  {"pgpgin",              &vm_pgpgin},     // important
  {"pgpgout",             &vm_pgpgout},     // important
  {"pgrefill",            &vm_pgrefill},  // GONE (now separate dma,high,normal)
  {"pgrefill_dma",        &vm_pgrefill_dma},
  {"pgrefill_high",       &vm_pgrefill_high},
  {"pgrefill_normal",     &vm_pgrefill_normal},
  {"pgrotated",           &vm_pgrotated},
  {"pgscan",              &vm_pgscan},  // GONE (now separate direct,kswapd and dma,high,normal)
  {"pgscan_direct_dma",   &vm_pgscan_direct_dma},
  {"pgscan_direct_high",  &vm_pgscan_direct_high},
  {"pgscan_direct_normal",&vm_pgscan_direct_normal},
  {"pgscan_kswapd_dma",   &vm_pgscan_kswapd_dma},
  {"pgscan_kswapd_high",  &vm_pgscan_kswapd_high},
  {"pgscan_kswapd_normal",&vm_pgscan_kswapd_normal},
  {"pgsteal",             &vm_pgsteal},  // GONE (now separate dma,high,normal)
  {"pgsteal_dma",         &vm_pgsteal_dma},
  {"pgsteal_high",        &vm_pgsteal_high},
  {"pgsteal_normal",      &vm_pgsteal_normal},
  {"pswpin",              &vm_pswpin},     // important
  {"pswpout",             &vm_pswpout},     // important
  {"slabs_scanned",       &vm_slabs_scanned},
  };
  const int vm_table_count = sizeof(vm_table)/sizeof(vm_table_struct);

#if __SIZEOF_LONG__ == 4
  unsigned long long slotll;
#endif

  vm_pgalloc = 0;
  vm_pgrefill = 0;
  vm_pgscan = 0;
  vm_pgsteal = 0;

  FILE_TO_BUF(VMINFO_FILE,vminfo_fd);

  head = buf;
  for(;;){
    tail = strchr(head, ' ');
    if(!tail) break;
    *tail = '\0';
    if(strlen(head) >= sizeof(namebuf)){
      head = tail+1;
      goto nextline;
    }
    strcpy(namebuf,head);
    found = bsearch(&findme, vm_table, vm_table_count,
        sizeof(vm_table_struct), compare_vm_table_structs
    );
    head = tail+1;
    if(!found) goto nextline;
#if __SIZEOF_LONG__ == 4
    // A 32 bit kernel would have already truncated the value, a 64 bit kernel
    // doesn't need to.  Truncate here to let 32 bit programs to continue to get
    // truncated values.  It's that or change the API for a larger data type.
    slotll = strtoull(head,&tail,10);
    *(found->slot) = (unsigned long)slotll;
#else
    *(found->slot) = strtoul(head,&tail,10);
#endif
nextline:

//if(found) fprintf(stderr,"%s=%d\n",found->name,*(found->slot));
//else      fprintf(stderr,"%s not found\n",findme.name);

    tail = strchr(head, '\n');
    if(!tail) break;
    head = tail+1;
  }
  if(!vm_pgalloc)
    vm_pgalloc  = vm_pgalloc_dma + vm_pgalloc_high + vm_pgalloc_normal;
  if(!vm_pgrefill)
    vm_pgrefill = vm_pgrefill_dma + vm_pgrefill_high + vm_pgrefill_normal;
  if(!vm_pgscan)
    vm_pgscan   = vm_pgscan_direct_dma + vm_pgscan_direct_high + vm_pgscan_direct_normal
                + vm_pgscan_kswapd_dma + vm_pgscan_kswapd_high + vm_pgscan_kswapd_normal;
  if(!vm_pgsteal)
    vm_pgsteal  = vm_pgsteal_dma + vm_pgsteal_high + vm_pgsteal_normal;
}


void cpuinfo (void) {
  // ought to count CPUs in /proc/stat instead of relying
  // on glibc, which foolishly tries to parse /proc/cpuinfo
  // note: that may have been the case but now /proc/stat
  //       is the default source.  parsing of /proc/cpuinfo
  //       only occurs if the open on /proc/stat fails
  //
  // SourceForge has an old Alpha running Linux 2.2.20 that
  // appears to have a non-SMP kernel on a 2-way SMP box.
  // _SC_NPROCESSORS_CONF returns 2, resulting in HZ=512
  // _SC_NPROCESSORS_ONLN returns 1, which should work OK

  smp_num_cpus = sysconf(_SC_NPROCESSORS_ONLN);
  if (smp_num_cpus<1)        /* SPARC glibc is buggy */
    smp_num_cpus=1;
}
