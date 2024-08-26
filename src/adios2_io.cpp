/*
 * Distributed under the OSI-approved Apache License, Version 2.0.  See
 * accompanying file Copyright.txt for details.
 *
 * hello-world.cpp : adios2 low-level API example to write and read a
 *                   std::string Variable with a greeting
 *
 *  Created on: Nov 14, 2019
 *      Author: William F Godoy godoywf@ornl.gov
 */
#include <Omega_h_timer.hpp>
#include <Omega_h_file.hpp>
#include <Omega_h_cmdline.hpp>
#include <iostream>
#include <stdexcept>

#include <adios2.h>
#if ADIOS2_USE_MPI
#include <mpi.h>
#endif

// to check the content of .bp, run /lore/seols/romulus-install/bin/bpls
// ex. ./bpls hello-world.bp
void writer(adios2::ADIOS &adios)
{ 
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    adios2::IO io = adios.DeclareIO("test-writer");

    const std::string greeting = "Hello World from ADIOS2";
    std::vector<float> myFloats;
    std::vector<int> myInts;
    for (int i=0; i<(rank+1); ++i)
    {
      myFloats.push_back(i*rank*2);
    }
    const std::size_t Nx = myFloats.size();
    myInts.push_back(Nx);

     adios2::Variable<float> bpFloats = io.DefineVariable<float>(
            "bpFloats", {size * Nx}, {rank * Nx}, {Nx}, adios2::ConstantDims);

     adios2::Variable<int> bpInts = io.DefineVariable<int>("bpInts", {size}, {rank},
		     {1}, adios2::ConstantDims);

    const std::string myString("Variable String from rank " + std::to_string(rank));

    for (int i=0; i<size; ++i)
    {
      if (rank == i)
        std::cout << "("<<rank<<") array size "<<Nx<<", string: "<<myString <<"\n";
      MPI_Barrier(MPI_COMM_WORLD);
    }

    adios2::Variable<std::string> bpString = io.DefineVariable<std::string>("bpString");
        (void)bpString;

    adios2::Variable<std::string> varGreeting = io.DefineVariable<std::string>("Greeting");

    std::string filename = "hello-world.bp";
    adios2::Engine writer = io.Open(filename, adios2::Mode::Write);
    writer.BeginStep();
    writer.Put(bpFloats, myFloats.data());
    writer.Put(bpInts, myInts.data());
    writer.Put(bpString, myString);
//    writer.Put(varGreeting, greeting);
    writer.EndStep();
    writer.Close();
    if (rank == 0)
        {
            std::cout << "Wrote file " << filename
                      << " to disk. It can now be read by running "
                         "adios2_install/bin/bpls.\n";
        }
}

void reader(adios2::ADIOS &adios)
{
    int rank, size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    adios2::IO io = adios.DeclareIO("hello-world-reader");
    adios2::Engine reader = io.Open("hello-world.bp", adios2::Mode::Read);

    reader.BeginStep();
    const std::map<std::string, adios2::Params> variables = io.AvailableVariables();
    for (const auto &variablePair : variables)
    {
      std::cout << "Name: " << variablePair.first;
      for (const auto &parameter : variablePair.second)
      {
        std::cout << "\t" << parameter.first << ": " << parameter.second << "\n";
      }
    }

    adios2::Variable<float> bpFloats = io.InquireVariable<float>("bpFloats");
    adios2::Variable<int> bpInts = io.InquireVariable<int>("bpInts");

    std::size_t Nx=1;

    if (bpInts) // means not found
    {
      std::vector<int> myInts;
      // read only the chunk corresponding to our rank
      bpInts.SetSelection({{rank}, {1}});

      reader.Get(bpInts, myInts, adios2::Mode::Sync);
      Nx=myInts[0];
    }

   if (bpFloats) // means found
  {
   std::vector<float> myFloats;

   // read only the chunk corresponding to our rank
   bpFloats.SetSelection({{Nx * rank}, {Nx}});
   reader.Get(bpFloats, myFloats, adios2::Mode::Sync);

   for (int i=0; i<size; ++i)
   { 
     if (rank == i)
     {
       std::cout << "array size "<<Nx<<"\n";
       for (const auto number : myFloats)
           std::cout <<"("<< rank<<") "<<number << " ";
     }
     MPI_Barrier(MPI_COMM_WORLD);
     std::cout << "\n";
   }
}

    adios2::Variable<std::string> varText = io.InquireVariable<std::string>("bpString");
    std::string text;
    reader.Get(varText, text);

//    adios2::Variable<std::string> varGreeting = io.InquireVariable<std::string>("Greeting");
//    std::string greeting;
//    reader.Get(varGreeting, greeting);

    for (int i=0; i<size; ++i)
    {
      if (rank == i)
        std::cout << "("<<rank<<") "<<text <<"\n";
      MPI_Barrier(MPI_COMM_WORLD);
    }
    reader.EndStep();
    reader.Close();
}

int main(int argc, char *argv[])
{
  auto lib = Omega_h::Library(&argc, &argv);
  auto t0 = Omega_h::now();
  Omega_h::Now t1;
  auto world = lib.world();
  auto rank = world->rank();
  auto size = world->size();

  std::vector<float> myFloats = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  std::vector<int> myInts = {0, -1, -2, -3, -4, -5, -6, -7, -8, -9};
  const std::size_t Nx = myFloats.size();

    const std::string myString("Hello Variable String from rank " + std::to_string(rank));

    try
    {
        adios2::ADIOS adios(MPI_COMM_WORLD);

	t0 = Omega_h::now();
        writer(adios);
	auto t1= Omega_h::now();
        if (!rank)
          std::cout << "write " << (t1 - t0)<<" (sec)\n";

        reader(adios);
	t0 = Omega_h::now();

        if (!rank) std::cout << "read " << (t0 - t1)<<" (sec)\n";
    }
    catch (std::exception &e)
    {
        std::cout << "ERROR: ADIOS2 exception: " << e.what() << "\n";
#if ADIOS2_USE_MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
    }
    return 0;
}
