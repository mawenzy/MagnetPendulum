// Main code changed to work with gnuplot

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>

#define RADIUS 0.5 // distance magnet - origin
#define h 0.02     // step size in calculation of differential equation
#define Z 0.25     // hight of the pendulum

/*
*   Function: getYForce
*   ------------------------------------
*   calculates the y-component of the force the magnets apply to the pendulum
*
*   (xalt , yalt):  position of the pendulum in the x-y plane
*   F:              Factor to change the strength of the force
*   N:              Amount of magnets, magnets are placed equally around a
*                       circle of radius 'RADIUS'
*
*   returns:        y - component of the force
*/
void getForce(double *force, double xalt, double yalt, double F, int N, double *x_val, double *y_val)
{
    int i;
    force[0] = 0;
    force[1] = 0;
    //iterate over all magnets
    for (i = 0; i < N; i++)
    {
        double magX = x_val[i]; //x coordinate of magnet
        double magY = y_val[i]; //y coordinate of magnet

        double dis = sqrt(Z * Z + (xalt - magX) * (xalt - magX) + (yalt - magY) * (yalt - magY));

        force[0] = force[0] + (xalt - magX) / (dis * dis * dis);
        force[1] = force[1] + (yalt - magY) / (dis * dis * dis);
    }
    force[0] = force[0] * F;
    force[1] = force[1] * F;
}

/*
*   Function: getColour
*   -------------------------------
*   evaluates the position of the pendulum, decides in which circle sector it is.
*
*   (x, y):         position in the x-y plane
*   N:              Amount of sectors, sectors are spaced with equal size around the origin
*
*   returns:        int corresponding to the magnet
*/
int getColour(double x, double y, int N)
{
    double angle = atan2(y, x);
    //special treatment of first sector, easy to calculate this way
    if (fabs(angle) < M_PI / N)
    {
        return 1;
    }
    //change to an angle between 0 and 2pi
    /*
        NOTE:
        To work with the angles in an easy way, you want 
            atan2( sin(t), cos(t)) = t for t in [0:2*pi],
        the next if statement fixes the result from the atan2
        function in a way making this equation true for all angles
    */
    if (angle < 0)
    {
        angle += 2 * M_PI;
    }
    int i;
    //iterator over the rest of the sectors
    for (i = 1; i < N; i++)
    {
        //test if the pendulum is in the sector
        if (angle >= 2 * M_PI * i / N - M_PI / N && angle <= 2 * M_PI * (i + 1) / N - M_PI / N)
        {
            return i + 1;
        }
    }
    return 0;
}

/*
*   Function: calcDiff
*   ----------------------------------------------
*   calculates the movement of the pendulum until it has stopped
*
*   (xalt, yalt):       starting position
*   F:                  coeffizient to the force
*   N:                  amount of magnets
*   Y:                  strength of drag
*   K:                  spring constant
*
*   returns:            int correspoding to the sector of the end position
*                           of the pendulum
*/
int calcDiff(double xalt, double yalt, double F, int N, double Y, double K, double *x_val, double *y_val)
{
    double vxalt = 0, vyalt = 0; // start values for speed

    double force[2] = {0};
    int amount = 0;
    int wasMoving = 0; //keeps track if the pendulum started moving, to start testing the stop condition (low cinetic energy) sooner

    //calculates the movement of the pendulum
    while (amount <= 5000) //default: 5000
    {
        //calculate acceleration
        getForce(force, xalt, yalt, F, N, x_val, y_val);
        double axalt = -K * xalt - Y * vxalt - force[0];
        double ayalt = -K * yalt - Y * vyalt - force[1];

        //get approximated position values using taylor expansion to the second term
        double xneu = xalt + h * vxalt + 0.5 * h * h * axalt;
        double yneu = yalt + h * vyalt + 0.5 * h * h * ayalt;

        //calculate new speed values
        double vxneu = vxalt + h * axalt;
        double vyneu = vyalt + h * ayalt;

        //shift position values
        xalt = xneu;
        yalt = yneu;
        vxalt = vxneu;
        vyalt = vyneu;

        //test if pendulum has started moving
        if (vxalt * vxalt + vyalt * vyalt > 0.1)
        {
            wasMoving = 1;
        }

        amount++;
        // at least 1000 iterations before testing the stop condition
        if (amount > 1000 || wasMoving == 1)
        {
            //energy small > pendulum is supposed to be in the final sector
            if (vxalt * vxalt + vyalt * vyalt < 0.01)
            {
                return getColour(xalt, yalt, N);
            }
        }
    }
    return getColour(xalt, yalt, N);
}

void gnuplot(FILE *gnuplot_pipe, double ZOOM, double STEP, double F, int N, double Y, double K, const char *outputFile)
{
    // commands to gnuplot
    fprintf(gnuplot_pipe, "set terminal png size 3000,3000\n");
    fprintf(gnuplot_pipe, "set output '%s'\n", outputFile);
    fprintf(gnuplot_pipe, "set margins 0,0,0,0\nunset xtics\nunset ytics\nunset border\nunset key\n");
    fprintf(gnuplot_pipe, "set palette model RGB maxcolors 10\n");
    fprintf(gnuplot_pipe, "set palette defined ( 0 'blue', 1 'yellow', 2 'red', 3 'green' )\n");
    fprintf(gnuplot_pipe, "unset colorbox\n");
    fprintf(gnuplot_pipe, "plot '-' matrix with image\n");

    //preparations
    double x_val[N];
    double y_val[N];
    int i;
    for (i = 0; i < N; i++)
    {
        x_val[i] = RADIUS * cos(2 * M_PI * i / N);
        y_val[i] = RADIUS * sin(2 * M_PI * i / N);
    }

    // go through a x-y grid
    for (double x = -ZOOM; x <= ZOOM + STEP; x += STEP)
    {
        for (double y = -ZOOM; y <= ZOOM + STEP; y += STEP)
        {
            //file output
            fprintf(gnuplot_pipe, "%d ", calcDiff(x, y, F, N, Y, K, x_val, y_val));
        }
        fprintf(gnuplot_pipe, "\n"); // newline after each row
    }

    // end gnuplot
    fprintf(gnuplot_pipe, "e\n");
}

int main(int argc, char **argv)
{
    double Y; //drag coefficient
    double K; //spring constant
    double F; //force of magnets
    int N;    //amount of magnets
    double ZOOM;
    double STEP = 0.005; // 0.005 -> 500x500
    char *outputFile = "magnet.png";

    //variables for reading in the input
    int qflag = 0, hflag = 0, oflag = 0, mflag = 0, zflag = 0, uflag = 0, fflag = 0, yflag = 0, kflag = 0;
    char *amountMagnets, *zoomString, *force, *yString, *kString;

    extern char *optarg;
    extern int optind;
    int c, err = 0;
    static char usage[] = "usage: magnets.exe [-qhu] [-o outputFile] [-f force] [-m magnets] [-z zoom] [-y drag] [-k spring_constant]\n";

    //o:outputFile m:magnets z:zoom f:force q_quiet h_help u_usage y:drag k:spring_constant o:outputfile
    while ((c = getopt(argc, argv, "o:m:z:f:qhuy:k:")) != -1)
    {
        switch (c)
        {
        case 'y':
            yflag = 1;
            yString = optarg;
            break;
        case 'k':
            kflag = 1;
            kString = optarg;
            break;
        case 'q':
            qflag = 1;
            break;
        case 'h':
            hflag = 1;
            break;
        case 'o':
            oflag = 1;
            outputFile = optarg;
            break;
        case 'm':
            mflag = 1;
            amountMagnets = optarg;
            break;
        case 'f':
            fflag = 1;
            force = optarg;
            break;
        case 'u':
            uflag = 1;
            break;
        case 'z':
            zflag = 1;
            zoomString = optarg;
        case '?':
            err = 1;
            break;
        }
    }

    //inpret input arguments
    //set amount of magnets
    if (mflag == 1)
    {
        N = atoi(amountMagnets);
    }
    else
    {
        N = 3;
    }
    //set zoom level
    if (zflag)
    {
        ZOOM = atof(zoomString);
    }
    else
    {
        ZOOM = 1.25;
    }
    //set force
    if (fflag)
    {
        F = atof(force);
    }
    else
    {
        F = 0.4;
    }
    //set drag
    if (yflag == 1)
    {
        Y = atof(yString);
    }
    else
    {
        Y = 0.2;
    }
    //set federkonstante
    if (kflag == 1)
    {
        K = atof(kString);
    }
    else
    {
        K = 0.3;
    }
    //print usage
    if (err == 1 || uflag == 1 || hflag == 1)
    {
        printf(usage);
        return 0;
    }

    //terminate program for false parameters, one magnet only would produce a one color image
    if (N == 1)
    {
        printf("really?");
        return 0;
    }

    //approx amount of pixels
    int PIXELS = (int)(2 * ZOOM / STEP) * (2 * ZOOM / STEP);

    //print info:
    if (qflag == 0)
    {
        printf("Amount of data points to calculate: %d\n", PIXELS);
        printf("Amount of magnets: %d, Force: %lg\n", N, F);
        printf("ZOOM: %lg\n", ZOOM);
        printf("Using gnuplot pipe\n");
        printf("Output file: %s\n\n", outputFile);
    }
    
    FILE *gnuplot_pipe = popen("gnuplot", "w");

    gnuplot(gnuplot_pipe, ZOOM, STEP, F, N, Y, K, outputFile);

    int close = fclose(gnuplot_pipe);
    if (close != 0)
    {
        printf("Error in closing the file");
        return 0;
    }
    return 0;
}