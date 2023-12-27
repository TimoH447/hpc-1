#include "mpi.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

double lapl(double x1, double x2,double x3,double x4,double x5)
{
    double a=-x1-x2-x3-x4+4.0*x5;
    return a;
}
double f(double x,double y){
    return 0;
}
double g(double x,double y){
    if (y==0){
        return 1-x;
    }
    else if (y==1){
        return x;
    }
    else if (x==0){
        return 1-y;
    }
    else {
        return y;
    }
}
double solution_u(double x, double y){
    return 2*x*y-x-y+1;
}

int main(int argc, char *argv[])
{
    //Für die Zeitnahme
    double t1,t2,t3,t4,t5;
    int i,j,k,size,rank;

    double w = 1.0;
    //Variablen für Abbruchkriterium
    double epsilon = 0.01; // Toleranz
    double squaredEuclideanError; 
    int loop_counter = 0;
    int MAX_LOOPS = 100;
    
    int print=0;
    //MPI Init
    MPI_Init(&argc, &argv);
    
    MPI_Request request;
    MPI_Status status;
    //Startzeit
    t1=MPI_Wtime();
    
    
    
    //Prozessnummer und Anzahl Prozesse holen
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int N=size;
    if (argc>1)
        N=atoi(argv[1]);
    
    if (argc>2)
        print=atoi(argv[2]);
    
    
    //Blockgrößen
    int n = N*N;
    int P = (int) sqrt(size);
    int m=N/size;
    int mP = N*P/size;
    printf("Blocksize: %d",mP);
    double h = (double) 1/ (N+1);
    
    double *x,*column;
    
    //y hat die Größe des Gebiets. x hat einen Layer drumrum, wo Nullen für den Rand und die kommunizierten Werte fürs Interface rein kommen.
    x=(double *)malloc((mP +2)*(mP+2)*sizeof(double));
    column = (double *)malloc((mP)*sizeof(double));
    
    //lokales x hat hier mal Nullen... kann man aber auch anders befüllen...
    for (i=1;i<mP+1;i++)
    {
        for (j=1;j<mP+1;j++)
        {
            x[i*(mP+2)+j]=0;
        }
    }
    
    // Rand befüllen mit nullen

    // Untere Prozessoren 1 bis P haben unten am Rand null einträge
    if (rank<P){
    //unten...
        for (i=0;i<mP+2;i++)
        {
            x[i]=g((rank*mP+i)*h,0.0);
        }
    }

    // Obere Prozessoren P^2 - P bis P^2 haben am oberen Rand null einträge
    if (rank>=size-P)
    {
        for (i=0;i<mP+2;i++)
        {
            x[(mP+1)*(mP+2)+i]=g((rank*mP+i)*h,1.0);
        }
    }

    // Linke Prozessoren Blöcke kP + 1 haben am linken Rand null einträge
    if ((rank+1)%P==1){
        for (i=0;i<mP+2;i++)
        {
            x[i*(mP+2)]=g(0,(rank/P)*mP+i);
        }
    }

    // Rechte Prozessoren Blöcke kP haben am rechten Rand null einträge
    if ((rank+1)%P==0){
        for (i=0;i<mP+2;i++)
        {
            x[i*(mP+2)+mP+1]=g(1,(((rank+1)/P)-1)*mP+i);
        }
    }

    // Tauschen der schwarzen Knoten am Partitionrand
    // (Wir tauschen alle Knoten, die paar mehr machen den Braten auch nicht fett)

    //Senden der lokalen "oberen Einträge im Streifen" von x nach oben, in den unteren Ghost Layer von x des Prozessors "darüber"
    // Prozessor soll keinen Gitterblock der obersten oder untersten Reihe sein
    if ((rank>=P) && (rank<size-P))
    {
        //Sende nach oben
        MPI_Isend(&x[mP*(mP+2)+1],mP, MPI_DOUBLE,rank+P, 41, MPI_COMM_WORLD, &request);
        
        //von unten empfangen
        MPI_Irecv(&x[(0)*(mP+2)+1], mP, MPI_DOUBLE, rank-P,41, MPI_COMM_WORLD, &request);
    }
    else if ((rank<P) && (size>P))
    {
        //Sende nach oben
        MPI_Isend(&x[mP*(mP+2)+1],mP, MPI_DOUBLE,rank+P, 41, MPI_COMM_WORLD, &request);
        
    }
    else if ((rank>=size-P) && (size>P))
    {
        //von unten empfangen
        MPI_Irecv(&x[(0)*(mP+2)+1], mP, MPI_DOUBLE, rank-P,41, MPI_COMM_WORLD, &request);
    }
    
    if (size>P)
        MPI_Wait(&request, &status);
    
    //Senden der lokalen "unteren Einträge im Streifen" von x nach unten, in den oberen Ghost Layer von x des Prozessors "darunter"
    // Prozessor soll keinen Gitterblock der obersten oder untersten Reihe sein
    if ((rank>=P) && (rank<size-P))
    {
        //Sende nach unten
        MPI_Isend(&x[(1)*(mP+2)+1],mP, MPI_DOUBLE,rank-P, 41, MPI_COMM_WORLD, &request);
        
        //von oben empfangen
        MPI_Irecv(&x[(mP+1)*(mP+2)+1], mP, MPI_DOUBLE, rank+P,41, MPI_COMM_WORLD, &request);
    }
    else if ((rank>=size-P) && (size>P))
    {
        //Sende nach unten
        MPI_Isend(&x[(1)*(mP+2)+1],mP, MPI_DOUBLE,rank-P, 41, MPI_COMM_WORLD, &request);
        
    }
    else if ((rank<P) && (size>P))
    {
        //von oben empfangen
        MPI_Irecv(&x[(mP+1)*(mP+2)+1], mP, MPI_DOUBLE, rank+P,41, MPI_COMM_WORLD, &request);
    }
    
    if (size>P)
        MPI_Wait(&request, &status);
    
    // Senden der lokalen "rechten Einträge" von x nach rechts in den rechten Ghost Layer von dem Prozessor mit dem Block rechts daneben
    if (((rank+1)%P>1)&&((rank+1)%P>0)){
        // Sende nach rechts
        // befüllen des column array (helper array for sending)
        for (i=1;i<mP+1;i++){
            column[i]=x[i*(mP+2)+mP];
        }
        MPI_Isend(&column,mP, MPI_DOUBLE,rank+1, 41, MPI_COMM_WORLD, &request);

        // Empfange von links
        MPI_Irecv(&column, mP, MPI_DOUBLE, rank-1,41, MPI_COMM_WORLD, &request);
        for (i=1;i<mP+1;i++){
            x[i*(mP+2)]=column[i];
        }

    }
    else if (((rank +1)%P==1)&&(size>P)){
        // Sende nach rechts
        for (i=1;i<mP+1;i++){
            column[i]=x[i*(mP+2)+mP];
        }
        MPI_Isend(&column,mP, MPI_DOUBLE,rank+1, 41, MPI_COMM_WORLD, &request);
    }
    else if (((rank+1)%P==0)&&(size>P)){
        // Empfange von links
        MPI_Irecv(&column, mP, MPI_DOUBLE, rank-1,41, MPI_COMM_WORLD, &request);
        for (i=1;i<mP+1;i++){
            x[i*(mP+2)]=column[i];
        }
    }
    if (size>P)
        MPI_Wait(&request, &status);

    // Und das ganze noch für die linken Einträge. Es gibt bestimmt eine bessere Variante als alle vier Fälle einzeln zu machen
    if (((rank+1)%P>1)&&((rank+1)%P>0)){
        // Sende nach links
        // befüllen des column array (helper array for sending)
        for (i=1;i<mP+1;i++){
            column[i]=x[i*(mP+2)+1];
        }
        MPI_Isend(&column,mP, MPI_DOUBLE,rank-1, 41, MPI_COMM_WORLD, &request);


        // Empfange von rechts
        MPI_Irecv(&column, mP, MPI_DOUBLE, rank+1,41, MPI_COMM_WORLD, &request);
        for (i=1;i<mP+1;i++){
            x[i*(mP+2)+mP+1]=column[i];
        }

    }
    else if (((rank+1)%P==0)&&(size>P)){
        // Sende nach links
        // befüllen des column array (helper array for sending)
        for (i=1;i<mP+1;i++){
            column[i]=x[i*(mP+2)+1];
        }
        MPI_Isend(&column,mP, MPI_DOUBLE,rank-1, 41, MPI_COMM_WORLD, &request);
    }
    else if (((rank +1)%P==1)&&(size>P)){
        // Empfange von rechts
        MPI_Irecv(&column, mP, MPI_DOUBLE, rank+1,41, MPI_COMM_WORLD, &request);
        for (i=1;i<mP+1;i++){
            x[i*(mP+2)+mP+1]=column[i];
        }
    }
    if (size>P){
        MPI_Wait(&request, &status);
    }
    

    // repeat until Abbruchbedingung
    // meine erste do while schleife in jeglicher Programmiersprache
    do {
    
        // Update red and then black (0=red, 1=black)
        //Blockweise Multiplikation.
        for (int color=0;color<2;color++){

            for (i=0;i<mP;i++)
            {
                for (j=0;j<mP;j++)
                {
                    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                    // Hier Multiplikation implementieren
                    // Nutzen Sie die funktion lapl (oben im Code)
                    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                    // only update colored knots
                    if ((i+j)%2==color){
                        x[(i+1)*(n+2)+j+1] = lapl(x[(i+1)*(n+2)+(j+1)-1],x[(i+1)*(n+2)+j+1+1],x[i*(n+2)+(j+1)],x[(i+2)*(n+2)+(j+1)],x[(i+1)*(n+2)+(j+1)]);
                    }
                }
            }

            /////////////////////////////////////
            // Tausche Rote/Schwarze Knoten am Partitionsrand (wir tauschen einfach immer alle knoten am Rand)
            //////////////////////////////////

            //Senden der lokalen "oberen Einträge im Streifen" von x nach oben, in den unteren Ghost Layer von x des Prozessors "darüber"
            // Prozessor soll keinen Gitterblock der obersten oder untersten Reihe sein
            if ((rank>=P) && (rank<size-P))
            {
                //Sende nach oben
                MPI_Isend(&x[mP*(mP+2)+1],mP, MPI_DOUBLE,rank+P, 41, MPI_COMM_WORLD, &request);
                
                //von unten empfangen
                MPI_Irecv(&x[(0)*(mP+2)+1], mP, MPI_DOUBLE, rank-P,41, MPI_COMM_WORLD, &request);
            }
            else if ((rank<P) && (size>P))
            {
                //Sende nach oben
                MPI_Isend(&x[mP*(mP+2)+1],mP, MPI_DOUBLE,rank+P, 41, MPI_COMM_WORLD, &request);
                
            }
            else if ((rank>=size-P) && (size>P))
            {
                //von unten empfangen
                MPI_Irecv(&x[(0)*(mP+2)+1], mP, MPI_DOUBLE, rank-P,41, MPI_COMM_WORLD, &request);
            }
            
            if (size>P)
                MPI_Wait(&request, &status);
            
            //Senden der lokalen "unteren Einträge im Streifen" von x nach unten, in den oberen Ghost Layer von x des Prozessors "darunter"
            // Prozessor soll keinen Gitterblock der obersten oder untersten Reihe sein
            if ((rank>=P) && (rank<size-P))
            {
                //Sende nach unten
                MPI_Isend(&x[(1)*(mP+2)+1],mP, MPI_DOUBLE,rank-P, 41, MPI_COMM_WORLD, &request);
                
                //von oben empfangen
                MPI_Irecv(&x[(mP+1)*(mP+2)+1], mP, MPI_DOUBLE, rank+P,41, MPI_COMM_WORLD, &request);
            }
            else if ((rank>=size-P) && (size>P))
            {
                //Sende nach unten
                MPI_Isend(&x[(1)*(mP+2)+1],mP, MPI_DOUBLE,rank-P, 41, MPI_COMM_WORLD, &request);
                
            }
            else if ((rank<P) && (size>P))
            {
                //von oben empfangen
                MPI_Irecv(&x[(mP+1)*(mP+2)+1], mP, MPI_DOUBLE, rank+P,41, MPI_COMM_WORLD, &request);
            }
            
            if (size>P)
                MPI_Wait(&request, &status);
            
            // Senden der lokalen "rechten Einträge" von x nach rechts in den rechten Ghost Layer von dem Prozessor mit dem Block rechts daneben
            if (((rank+1)%P>1)&&((rank+1)%P>0)){
                // Sende nach rechts
                // befüllen des column array (helper array for sending)
                for (i=1;i<mP+1;i++){
                    column[i]=x[i*(mP+2)+mP];
                }
                MPI_Isend(&column,mP, MPI_DOUBLE,rank+1, 41, MPI_COMM_WORLD, &request);

                // Empfange von links
                MPI_Irecv(&column, mP, MPI_DOUBLE, rank-1,41, MPI_COMM_WORLD, &request);
                for (i=1;i<mP+1;i++){
                    x[i*(mP+2)]=column[i];
                }

            }
            else if (((rank +1)%P==1)&&(size>P)){
                // Sende nach rechts
                for (i=1;i<mP+1;i++){
                    column[i]=x[i*(mP+2)+mP];
                }
                MPI_Isend(&column,mP, MPI_DOUBLE,rank+1, 41, MPI_COMM_WORLD, &request);
            }
            else if (((rank+1)%P==0)&&(size>P)){
                // Empfange von links
                MPI_Irecv(&column, mP, MPI_DOUBLE, rank-1,41, MPI_COMM_WORLD, &request);
                for (i=1;i<mP+1;i++){
                    x[i*(mP+2)]=column[i];
                }
            }
            if (size>P)
                MPI_Wait(&request, &status);

            // Und das ganze noch für die linken Einträge. Es gibt bestimmt eine bessere Variante als alle vier Fälle einzeln zu machen
            if (((rank+1)%P>1)&&((rank+1)%P>0)){
                // Sende nach links
                // befüllen des column array (helper array for sending)
                for (i=1;i<mP+1;i++){
                    column[i]=x[i*(mP+2)+1];
                }
                MPI_Isend(&column,mP, MPI_DOUBLE,rank-1, 41, MPI_COMM_WORLD, &request);


                // Empfange von rechts
                MPI_Irecv(&column, mP, MPI_DOUBLE, rank+1,41, MPI_COMM_WORLD, &request);
                for (i=1;i<mP+1;i++){
                    x[i*(mP+2)+mP+1]=column[i];
                }

            }
            else if (((rank+1)%P==0)&&(size>P)){
                // Sende nach links
                // befüllen des column array (helper array for sending)
                for (i=1;i<mP+1;i++){
                    column[i]=x[i*(mP+2)+1];
                }
                MPI_Isend(&column,mP, MPI_DOUBLE,rank-1, 41, MPI_COMM_WORLD, &request);
            }
            else if (((rank +1)%P==1)&&(size>P)){
                // Empfange von rechts
                MPI_Irecv(&column, mP, MPI_DOUBLE, rank+1,41, MPI_COMM_WORLD, &request);
                for (i=1;i<mP+1;i++){
                    x[i*(mP+2)+mP+1]=column[i];
                }
            }
            if (size>P){
                MPI_Wait(&request, &status);
            }
            // Knoten Tauschende
        }

        // Update Abbruchkriterium
        // compute local squared error norm
        loop_counter+=1;
        double local_squaredEuclideanError=0;
        for (i=0;i<mP;i++)
        {
            for (j=0;j<mP;j++)
            {
                local_squaredEuclideanError += (x[(i+1)*(n+2)+j+1]-solution_u((double)i*h,(double)j*h))*(x[(i+1)*(n+2)+j+1]-solution_u((double)i*h,(double)j*h)); 
            }
        }
        MPI_Allreduce(&local_squaredEuclideanError,&squaredEuclideanError,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    
    } while ((sqrt(squaredEuclideanError)>epsilon) || (loop_counter > MAX_LOOPS));


    t2=MPI_Wtime();

    
    
    if (rank==0)
        printf("\n\n Benoetigte Zeit: %fs\n\n",t2-t1);
        printf("Euklidische Norm des Fehlervektors: %f\n",squaredEuclideanError);

    free(x);
    free(column);
    MPI_Finalize();
    return 0;
}