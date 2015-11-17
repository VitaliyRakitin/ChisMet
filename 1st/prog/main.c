//
//  main.c
//  
//
//  Created by Виталий Ракитин on 25.09.15.
//
//

#include "task1.h"

int  main(void)
{   int errflag = 0;
    double local_7[CHECK_POSITIONS * VAL_NUM];
    double local_9[CHECK_POSITIONS * VAL_NUM];
    double local_11[CHECK_POSITIONS * VAL_NUM];
    double check = 0;

    double global_error = 0;

    if (errflag == 0)
        errflag = runge(&global_error, 1e-7, "data_-7.csv", local_7);
    printf("Global error for -7: %e\n", global_error);
    if (errflag == 0)
        errflag = runge(&global_error, 1e-9, "data_-9.csv", local_9);
    printf("Global error for -9: %e\n", global_error);
    
    if (errflag == 0)
        errflag = runge(&global_error, 1e-11,"data_-11.csv", local_11);
    printf("Global error for -11: %e\n", global_error);

    printf("\n Local errors for:\n");
    for (int i=0;i<CHECK_POSITIONS*VAL_NUM;i++){
        if ((local_7[i] != 0)&&(local_9[i] != 0)&&(local_11[i] != 0)){
            printf("%d",i);
            check = (local_7[i] - local_9[i])/(local_9[i] - local_11[i]);
            printf("local_error = %e\n", check);
        }
    }

    if (errflag == -1) printf("File error!\n");
    if (errflag == -2) printf("Default too small step!\n");
    if (errflag == -3) printf("Memory error!\n");

    return 0;
}

