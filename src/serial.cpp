#include "serial.hpp"
#include "display.hpp"

/**
 * @brief Entry point for serial implementation
 */
void serial_main(void) {
    fprintf(stderr, "Hello serial\n");

    while (1) {
        display_render();
    }
}
