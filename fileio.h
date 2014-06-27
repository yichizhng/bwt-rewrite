#ifndef _FILEIO_H
#define _FILEIO_H

void write_index(const fm_index *fmi, FILE *f);

fm_index *read_index(FILE *f);

#endif /* _FILEIO_H */
