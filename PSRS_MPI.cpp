using namespace std;
#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <algorithm>
#include <string.h>
#include <vector>
#include <queue>
#include <math.h>
#include "mpi.h"

// compare_ints() used by cstdlib::qsort() to compare two integers
int compara(const void *i, const void *j)
{
	int int1 = *reinterpret_cast<const int *>(i);
	int int2 = *reinterpret_cast<const int *>(j);
	if (int1<int2) return -1;
	if (int1>int2) return 1;
	return 0;
}

// check array ordinato
string checkSort(int array[], int arrayStart, int arrayLength)
{
	for(int i=arrayStart; i<arrayStart+arrayLength-1;i++)
	{
		if (array[i]>array[i+1])
			return "non è ";
	}

	return " è ";
}

string checkSort(int array[], int arrayLength)
{
	return checkSort(array,0,arrayLength);
}


//merging
// This structure is what is placed in the priority queue: an index to the
//		appropriate array,  an index to an element in the appropriate array, and
//		the value of stored at the index of the element in the appropriate array
struct m_data
{
	int st_index;
	int index;
	int st_value;

	m_data(int st=0, int id=0, int stv = 0):st_index(st),index(id),st_value(stv){}

};


// comparison operator
bool operator<( const m_data & uno, const m_data & due)
{
	return uno.st_value > due.st_value;
}

int multimerge(int * starts[], const int len[], const int Number,int newArray[], const int newArray_Length)
{
 	// Create priority queue.  There will be at most one item in the priority queue
 	// for each of the Number lists.
 	priority_queue< m_data> priorit;

 	// Examine each of the Number start[] lists, place the first location into
	// the priority 	queue if the list is not empty
 	for(int i=0; i<Number;i++)
 	{
		if (len[i]>0)
		{
			priorit.push(m_data(i,0,starts[i][0]));
		}
	}
	// As long as priorities is not empty, pull off the top member (the smallest
	//value from list i), push it into the newArray, and place the next element from
	// list i in the priority queue
	int newArray_index = 0;  // index into the merged array
	while (!priorit.empty() && (newArray_index<newArray_Length))
	{
		// grab the smallest element, and remove it from the priority queue
		m_data array = priorit.top();
		priorit.pop();

		// insert this smallest element into the merged array
		newArray[newArray_index++] = starts[array.st_index][array.index];

		// if start[array.st_index] is not empty, place the next member into priority
		if ( len[array.st_index]>(array.index+1))
		{
			priorit.push(m_data(array.st_index, array.index+1,
								starts[array.st_index][array.index+1]));
		}
}

// return the logical size of the merged array
return newArray_index;
}


//****  MAIN

// La funzione main viene eseguita da ogni singolo processore
int main( int argc, char* argv[])
{
	// identificativo del processore (rank), numero totali di processori
    int IDproc, num_procs;

	//variabile usata per il tempo dal processore principale (#0)
    double startwtime = 0.0, endwtime;


	// *******************************************
	//
	//
	// FASE 0 : Inizializzazione

	// facciamo partire MPI e passiamo la linea di comando a MPI
    MPI::Init(argc, argv);

	// raccogliamo informazioni riguardo i
	// il numero di processori e l'identificativo di questo processore
    num_procs = MPI::COMM_WORLD.Get_size();
    IDproc = MPI::COMM_WORLD.Get_rank();

	// look through arguments for:
	//		 -DS nnnn to set myData_Size
	//		 -SR nnnn to set seed_rand
	int seed_rand = 1000;
	int myData_Size = 4000;
	for(int i=0;i<argc;i++)
	{
		// search for special arguments
		if (strcmp(argv[i],"-DS")==0)
		{
			myData_Size = atoi(argv[i+1]); i++;
		}
		else if (strcmp(argv[i],"-SR")==0)
		{
			seed_rand = atoi(argv[i+1]); i++;
		}
	}

    // array per immagazzinare i dati iniziali non ordinati
	// myData_Lengths[] e myData_Starts[] usati con Scatterv() per distribuire
	// myData[] agli altri processori
	int myData[myData_Size];
	int myData_Lengths[num_procs];
	int myData_Starts[num_procs];

	// buffer di comunicazione usato per la determinazione dei valori dei pivot
	int pivot_buff[num_procs*num_procs];
	int pivot_buffSize;


	// Calcoliamo la lunghezza individuale di mydata[] da esser distribuita
	// a processore numproc. L'ultimo processore avrà alcuni "extras"
	for(int i=0;i<num_procs;i++)
	{
		myData_Lengths[i] = myData_Size/num_procs;
		myData_Starts[i]= i*myData_Size/num_procs;
	}
	myData_Lengths[num_procs-1]+=(myData_Size%num_procs);


	// Il nodo principale inizializza i dati da testare, e fa anche partire il timer
    if (IDproc == 0)
    {
		// seme random
		srandom(seed_rand);
		for(int index=0; index<myData_Size; index++)
		{
			myData[index] = random()% 900;
		}

		startwtime = MPI::Wtime();
    }


	// *******************************************
	//
	// FASE I:  Distribuiamo i dati, ordinamento locale e raccolta di campioni

	// I dati sono distribuiti a tutti i processori dal processore principale (#0)
	// Il processore principale esegue "in place"
	if (IDproc==0)
	{
		MPI::COMM_WORLD.Scatterv(myData,myData_Lengths,myData_Starts,MPI::INT,MPI_IN_PLACE,myData_Lengths[IDproc],MPI::INT,0);
	}
	else
	{
		MPI::COMM_WORLD.Scatterv(myData,myData_Lengths,myData_Starts,MPI::INT,myData,myData_Lengths[IDproc],MPI::INT,0);
	}

	// Tutti i processori ordinano il proprio pezzo di dati usando quicksort
	qsort(myData,myData_Lengths[IDproc], sizeof(int), compara);


	//Tutti i processori raccolgono campioni regolari dalla lista ordinata considerando un offset pari all'indice di  myData[]
	for(int index=0;index<num_procs;index++)
	{
		pivot_buff[index]= myData[index*myData_Lengths[IDproc]/num_procs];
	}


	// *******************************************
	//
	// FASE II:  Raccolta e fusione dei campioni, e trasmissione di p-1 pivot


	// Il processore principale raccoglie tutti i candidati pivot dai processori, il processore principale ha in dati "in place"

	if (IDproc==0)
	{
		MPI::COMM_WORLD.Gather(MPI_IN_PLACE,num_procs,MPI::INT,
			pivot_buff,num_procs,MPI::INT,0);
	}
	else
	{
		MPI::COMM_WORLD.Gather(pivot_buff,num_procs,MPI::INT,
			pivot_buff,num_procs,MPI::INT,0);
	}

	// Il processore principale unisce le liste insieme e selezione
	// il valori del pivot normale che devono essere trasmessi
	if (IDproc == 0)
	{

		//unisce la lista del processore num_procs in una sola
		int *starts[num_procs];  //array delle liste
		int lengths[num_procs];  //array della lunghezza delle liste
		for(int i=0;i<num_procs;i++)
		{
			starts[i]=&pivot_buff[i*num_procs];
			lengths[i]=num_procs;
		}
		int tempbuffer[num_procs*num_procs];  // lista "mergiata"
		multimerge(starts,lengths,num_procs,tempbuffer,num_procs*num_procs);

		//seleziona regolarmente numprocs-1 pivot da trasmettere
		// come partizione dei valori di pivot a myData
		for(int i=0; i<num_procs-1; i++)
		{
			pivot_buff[i] = tempbuffer[(i+1)*num_procs];
		}
	}

	// Root processor (#0) broadcasts the partition values
	MPI::COMM_WORLD.Bcast(pivot_buff,num_procs-1,MPI::INT,0);


	// *******************************************
	//
	// FASE III: Partizione dei dati locali

	// Tutti i processori partizionano la loro parte di dati basandosi
	// sui valori del pivot memorizzati in pivot_buff[]

	// Informazioni di partizionamento per myData[]:
	// 				l'indice dell'inizione della i-esaima classe è class_Start[i]
	//				la lunghezza della i-esima classe è classLenght[i]
	//				e i membri della i-esima classe, myData[j], hanno la seguente proprietà
	//				  pivot_buff[i-1]<= myData[j] < pivot_buff[i]

	int class_Start[num_procs];
	int class_Length[num_procs];


	// è necessario  che ogni processore partizioni la propria lista usando i valori di pivot_buffs
	int dataindex=0;
	for(int classindex=0; classindex<num_procs-1; classindex++)
	{
		class_Start[classindex] = dataindex;
		class_Length[classindex]=0;

		// purché dataindex si riferisca ai dati nella classe corrente
		while((dataindex< myData_Lengths[IDproc])
			&& (myData[dataindex]<=pivot_buff[classindex]))
		{
			class_Length[classindex]++;
			dataindex++;
		}
	}
	// imposta start and lenght per l'ultima classe
	class_Start[num_procs-1] = dataindex;
	class_Length[num_procs-1] = myData_Lengths[IDproc] - dataindex;


	// *******************************************
	//
	// FASE IV: Tutte le i-esima classi sono raccolte dal processore i
	int recv_buff[myData_Size];    //buffer per tenere in memoria tutti i membri della classe i
	int recv_Lengths[num_procs];     // numero del processore, lunghezza di ogni numerodiprocessore-iesima classe
	int recv_Starts[num_procs];      //indice che ci dice dove iniziare a memorizzare


	//si prendono tutte le liste orindare e tutti i valori ordinati del iprocessore-esima classe
	for(int iprocessor=0; iprocessor<num_procs; iprocessor++)
	{
		// Each processor, iprocessor gathers up the numproc lengths of the sorted
		// values in the iprocessor class
		MPI::COMM_WORLD.Gather(&class_Length[iprocessor], 1, MPI::INT,recv_Lengths,1,MPI::INT,iprocessor);


		// From these lengths the myid^th class starts are computed on
		// processor IDproc
		if (IDproc == iprocessor)
		{
			recv_Starts[0]=0;
			for(int i=1;i<num_procs; i++)
			{
				recv_Starts[i] = recv_Starts[i-1]+recv_Lengths[i-1];
			}
		}

		// each iprocessor gathers up all the members of the iprocessor^th
		// classes from the other nodes
		MPI::COMM_WORLD.Gatherv(&myData[class_Start[iprocessor]],class_Length[iprocessor],MPI::INT,
			recv_buff,recv_Lengths,recv_Starts,MPI::INT,iprocessor);
	}


	// multimerge these numproc lists on each processor
	int *m_Starts[num_procs]; // array of list starts
	for(int i=0;i<num_procs;i++)
	{
		m_Starts[i]=recv_buff+recv_Starts[i];
	}
	multimerge(m_Starts,recv_Lengths,num_procs,myData,myData_Size);

	int mysend_Length = recv_Starts[num_procs-1] + recv_Lengths[num_procs-1];

	// *******************************************
	//
	// PHASE VI:  Root processor collects all the data


	int send_Lengths[num_procs]; // lengths of consolidated classes
	int send_Starts[num_procs];  // starting points of classes
	// Root processor gathers up the lengths of all the data to be gathered
	MPI::COMM_WORLD.Gather(&mysend_Length,1,MPI::INT,send_Lengths,1,MPI::INT,0);

	// The root processor compute starts from lengths of classes to gather
	if (IDproc == 0)
	{
		send_Starts[0]=0;
		for(int i=1; i<num_procs; i++)
		{
			send_Starts[i] = send_Starts[i-1]+send_Lengths[i-1];
		}
	}

	// Now we let processor #0 gather the pieces and glue them together in
	// the right order
	int arraySorted[myData_Size];
	MPI::COMM_WORLD.Gatherv(myData,mysend_Length,MPI::INT,arraySorted,send_Lengths,send_Starts,MPI::INT,0);

	//tempo sequenziale caso medio O(n^2) con n numero degli elementi
    double tseq = myData_Size * double(log10(myData_Size)) ;
    double speed_up = tseq / -(endwtime-startwtime) ;
    double efficienza = speed_up / num_procs;

    if (IDproc == 0)
	{
		endwtime = MPI::Wtime();
        cout << "\nTempo di esecuzione totale (in secondi) = "
		     <<  setprecision(4) << endwtime-startwtime << endl;

		cout << "\nESITO TEST: il data set " << checkSort(arraySorted,myData_Size) << " ordinato!" << endl;
/*
		cout << "\nPRESTAZIONI:\n"
		     <<  "\tSPEED-UP = " << speed_up << "\n\tEFFICIENZA = " << efficienza <<"\n"<< endl;

*/
	}

    // shutdown MPI on the processor
    MPI::Finalize();
    return 0;
}
