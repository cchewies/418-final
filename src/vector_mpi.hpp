#pragma once

#include "stars.hpp"

void mpi_send_stars(int receiver, std::vector<Star> stars);
std::vector<Star> mpi_recv_stars(int sender);
