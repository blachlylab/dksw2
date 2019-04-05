module dksw2.kalloc;

import core.stdc.stdlib;

// #define km_size(x) (*(((size_t*)(x))-1) * sizeof(size_t))

extern(C):

void *kmalloc(void *km, size_t size);
void *krealloc(void *km, void *ptr, size_t size);
void *kcalloc(void *km, size_t count, size_t size);
void kfree(void *km, void *ptr);

void *km_init();
void km_destroy(void *km);

void km_stat(const(void) *km); // TODO: return numbers instead of print to stderr

