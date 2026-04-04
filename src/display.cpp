#include "display.hpp"
#include <SDL2/SDL.h>
#include <iostream>

// SDL objects
static SDL_Window* sdlwindow = nullptr;
static SDL_Renderer* sdlrenderer = nullptr;

/**
 * @brief Handle keyboard input
 * 
 * @param event SDL2 event
 */
static void handle_key(const SDL_Event& event) {
    if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_q) {
        fprintf(stderr, "Quit");
        exit(0);
    }
}

/**
 * @brief Render a frame
 * 
 */
#include <cmath> // for log

void display_render(std::vector<Star> &stars) {

    // Handle key events
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        handle_key(event);
    }

    // Clear screen
    SDL_SetRenderDrawColor(sdlrenderer, 25, 50, 75, 255); // dark bluish background
    SDL_RenderClear(sdlrenderer);

    // Draw stars with size proportional to log(mass)
    SDL_SetRenderDrawColor(sdlrenderer, 255, 255, 255, 255); // white color

    for (auto& s : stars) {
        // Ensure mass is positive to avoid log(0)
        float mass = std::max(0.1f, s.mass);

        // Size based on log of mass, scaled
        int size = std::max(1, static_cast<int>(std::log(mass + 1.0f) * 3.0f)); 

        SDL_Rect rect;
        rect.x = static_cast<int>(s.x) - size / 2;
        rect.y = static_cast<int>(s.y) - size / 2;
        rect.w = size;
        rect.h = size;

        SDL_RenderFillRect(sdlrenderer, &rect);
    }

    SDL_RenderPresent(sdlrenderer);
}

/**
 * @brief Initialize display
 */
void display_init(void) {
    SDL_Init(SDL_INIT_VIDEO);

    sdlwindow = SDL_CreateWindow("barnes hut",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        WINDOW_WIDTH, WINDOW_HEIGHT, 0);

    sdlrenderer = SDL_CreateRenderer(sdlwindow, -1, SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
}
