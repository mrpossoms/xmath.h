#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <inttypes.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>

#define TEST int main (int argc, const char* argv[])

#define EPHEM_EPISLON 2e-9 // Allow for mantissa error

#define UT_CHECK_DOUBLE_TOL(val, expected, msg, tol)                          \
  if (fabs (val - expected) >= tol)                                           \
  {                                                                           \
    printf ("delta %s = %.5e\tval = %.16f\txpt = %.16f\n", msg,               \
            fabs (val - expected), val, expected);                            \
    assert (false);                                                           \
  }                                                                           \
  do                                                                          \
  {                                                                           \
  } while (0)

#define UT_CHECK_STATE(state, expected_state, num, tolerance)                 \
  UT_CHECK_DOUBLE_TOL (state.position.v[0], expected_state.position.v[0],         \
                       "ephem " #num " position[0]", tolerance);              \
  UT_CHECK_DOUBLE_TOL (state.position.v[1], expected_state.position.v[1],         \
                       "ephem " #num " position[1]", tolerance);              \
  UT_CHECK_DOUBLE_TOL (state.position.v[2], expected_state.position.v[2],         \
                       "ephem " #num " position[2]", tolerance);              \
  UT_CHECK_DOUBLE_TOL (state.velocity.v[0], expected_state.velocity.v[0],         \
                       "ephem " #num " velocity[0]", tolerance);              \
  UT_CHECK_DOUBLE_TOL (state.velocity.v[1], expected_state.velocity.v[1],         \
                       "ephem " #num " velocity[1]", tolerance);              \
  UT_CHECK_DOUBLE_TOL (state.velocity.v[2], expected_state.velocity.v[2],         \
                       "ephem " #num " velocity[2]", tolerance);              \
  UT_CHECK_DOUBLE_TOL (state.epoch_jde_days, expected_state.epoch_jde_days,   \
                       "ephem " #num " dt", tolerance)

void*
my_get_file_buffer (const char* path, size_t* buf_size, uint64_t* ckhsum)
{
  static void* mapped_file_buf;
  static size_t last_size;

  if (mapped_file_buf)
  {
    munmap (mapped_file_buf, last_size);
    mapped_file_buf = NULL;
  }

  if (NULL == path)
  {
    return NULL;
  }
  if (NULL == buf_size)
  {
    return NULL;
  }
  if (NULL == ckhsum)
  {
    return NULL;
  }

  // open the file
  int fd = open (path, O_RDONLY);
  if (fd < 0)
  {
    return NULL;
  }

  // get the size
  *buf_size = lseek (fd, 0, SEEK_END);
  lseek (fd, 0, SEEK_SET);
  last_size = *buf_size;

  // map the file
  mapped_file_buf = mmap (NULL, last_size, PROT_READ, MAP_PRIVATE, fd, 0);

  return mapped_file_buf;
}