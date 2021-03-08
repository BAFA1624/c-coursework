#include <stdio.h>
#include <stdlib.h>
#include <string.h>

size_t strlcat(char* dest, const char* src, size_t size);

int main(int argc, char* argv[])
{
    int status;
    char* my_err_msg = "\nUnknown error\n";
    char* push_msg = "\nDefault msg\n";
    if (argc > 1) {
	status = system("git add ../23762.c ../CMakeLists.txt ../final-exercise2020-2021.pdf ../README.md ../Plots 3_b_1.txt 3_b_2.txt git_push.c inv_1.txt inv_2.txt plot_3c.py") / 256;
	if (status != 0) {
	    my_err_msg = "\nCommand \"git add ....\" failed.\n";
	    goto error;
	}
	// This is a test
	// 14
	// git commit -m \"Auto Push: \"
	// 29
	size_t sz = strlen(argv[1]) + 29;
	printf("\nsz = %ld\n", sz);
	push_msg = (char*)malloc(sz * sizeof(char));
	if (!push_msg) {
	    my_err_msg = "\nFailed mem alloc for push_msg\n";
	    goto error;
	}
	size_t x = strlcat(push_msg, "git commit -m \"Auto Push: ", sz);
	printf("\n--\n%s|\n--\n", push_msg);
	size_t y = strlcat(push_msg, argv[1], sz);
	printf("\n--\n%s|\n--\n", push_msg);
	size_t z = strlcat(push_msg, "\"", sz);

	printf("\n--\n%s|\n--\n", push_msg);
	printf("\nx = %ld, y = %ld, z = %ld\n", x, y, z);

	status = system(push_msg) / 256;
	if (status != 0) {
	    my_err_msg = "\nCommand \"git commit ...\" failed.\n";
	    goto error;
	}

	status = system("git push") / 256;
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
