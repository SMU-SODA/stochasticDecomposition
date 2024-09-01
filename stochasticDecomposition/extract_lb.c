#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINE_LENGTH 1024
#define TARGET_LINE "Lower bound estimate"

int main(int argc, char *argv[])
{
    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s <output_file>\n", argv[0]);
        return 1;
    }

    FILE *file = fopen(argv[1], "r");
    if (!file)
    {
        perror("Error opening file");
        return 1;
    }

    char line[MAX_LINE_LENGTH];
    while (fgets(line, sizeof(line), file))
    {
        if (strstr(line, TARGET_LINE) != NULL)
        {
            // Extract the value after the colon
            strtok(line, ":");
            char *lb_str = strtok(NULL, ":");
            double lb_estimate = atof(lb_str);
            printf("%.6lf\n", lb_estimate); // Print with 6 decimal places
            break;
        }
    }

    fclose(file);
    return 0;
}
