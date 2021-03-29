#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef HAVE_STRLCAT
/*
 * '_cups_strlcat()' - Safely concatenate two strings.
 */

size_t		     /* O - Length of string */
strlcat(char* dst,   /* O - Destination string */
    const char* src, /* I - Source string */
    size_t size)     /* I - Size of destination string buffer */
{
    size_t srclen; /* Length of source string */
    size_t dstlen; /* Length of destination string */

    /*
  * Figure out how much room is left...
  */

    dstlen = strlen(dst);
    size -= dstlen + 1;

    if (!size)
	return (dstlen); /* No room, return immediately... */

    /*
  * Figure out how much room is needed...
  */

    srclen = strlen(src);

    /*
  * Copy the appropriate amount...
  */

    if (srclen > size)
	srclen = size;

    memcpy(dst + dstlen, src, srclen);
    dst[dstlen + srclen] = '\0';

    return (dstlen + srclen);
}
#endif /* !HAVE_STRLCAT */

#ifndef HAVE_STRLCPY
/*
 * '_cups_strlcpy()' - Safely copy two strings.
 */

size_t		     /* O - Length of string */
strlcpy(char* dst,   /* O - Destination string */
    const char* src, /* I - Source string */
    size_t size)     /* I - Size of destination string buffer */
{
    size_t srclen; /* Length of source string */

    /*
  * Figure out how much room is needed...
  */

    size--;

    srclen = strlen(src);

    /*
  * Copy the appropriate amount...
  */

    if (srclen > size)
	srclen = size;

    memcpy(dst, src, srclen);
    dst[srclen] = '\0';

    return (srclen);
}
#endif /* !HAVE_STRLCPY */

int main(int argc, char* argv[])
{
    int status;
    char* my_err_msg = "\nUnknown error\n";
    char* push_msg = "\nDefault msg\n";
    if (argc > 1) {
	status = system("git add ../23762.c ../CMakeLists.txt ../final-exercise2020-2021.pdf ../README.md ../Plots git_push.c h1.txt h2.txt inv_1.txt inv_2.txt inv_3.txt plot_3c.py plot_3h.py plot_3m.py") / 256;
	if (status != 0) {
	    my_err_msg = "\nCommand \"git add ...\" failed.\n";
	    goto error;
	}

	size_t sz = strlen(argv[1]) + 29;

	push_msg = (char*)malloc(sz * sizeof(char));
	if (!push_msg) {
	    my_err_msg = "\nFailed mem alloc for push_msg\n";
	    goto error;
	}

	strlcat(push_msg, "git commit -m \"Auto Push: ", sz);
	strlcat(push_msg, argv[1], sz);
	strlcat(push_msg, "\"", sz);

	status = system(push_msg) / 256;
	if (status != 0) {
	    my_err_msg = "\nCommand \"git commit ...\" failed.\n";
	    goto error;
	}

	status = system("git push -f") / 256;
	if (status != 0) {
	    my_err_msg = "Command \"git push\" failed.\n";
	    goto error;
	}

    } else {
	my_err_msg = "\nError likely due too missing cmd line args.\n";
	goto error;
    }

    goto cleanup;

error:
    fprintf(stderr, "\nAn error occurred.\nEnsure that command is of the form './git_push \"<your message here>\"'.\nErrorMsg: %s", my_err_msg);
    exit(-1);

cleanup:
    printf("\nGit push successful! Cleaning up...\n");
    free(push_msg);
    printf("\nFinished!\n");
}
