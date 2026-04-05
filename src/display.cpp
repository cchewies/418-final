/**
 * @file display.cpp
 * @author Zhuoyi Zou (zhuoyiz@andrew.cmu.edu)
 * 
 * Display functions for rendering output
 * Uses SDL2 and is kinda slow over AFS, probably good to not
 * render every single frame
 */

#include "display.hpp"
#include <SDL2/SDL.h>
#include <iostream>
#include <cmath>

// SDL objects
static SDL_Window* sdlwindow = nullptr;
static SDL_Renderer* sdlrenderer = nullptr;

static int frame_ctr = -1;

/**
 * @brief Poll for quit key
 * @return If quit command sent
 */
bool check_quit(void) {
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_q) {
            fprintf(stderr, "Quit");
            return true;
        }
    }
    return false;
}

/**
 * @brief Render a frame
 */
void display_render(std::vector<Star> &stars) {

#ifdef RENDER_ENABLED

    frame_ctr++;
    if (frame_ctr % RENDER_PERIOD != 0) {
        fprintf(stdout, "Render iteration %d, skipping\n", frame_ctr);
        return;
    }
    fprintf(stdout, "Render iteration %d, rendering\n", frame_ctr);

    // Clear screen
    SDL_SetRenderDrawColor(sdlrenderer, 25, 50, 75, 255); // dark bluish background
    SDL_RenderClear(sdlrenderer);

    SDL_SetRenderDrawColor(sdlrenderer, 255, 255, 255, 255); // white color

    for (auto& s : stars) {
        // Draw stars with size proportional to log(mass)
        float mass = std::max(0.0f, s.mass);
        int size = std::max(1, (int)(std::log(mass + 1) * 3)); 

        SDL_Rect rect;
        rect.x = static_cast<int>(s.x) - size / 2;
        rect.y = static_cast<int>(s.y) - size / 2;
        rect.w = size;
        rect.h = size;

        SDL_RenderFillRect(sdlrenderer, &rect);
    }

    SDL_RenderPresent(sdlrenderer);

#endif /* RENDER_ENABLED */

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

/**
 * @brief Cleanup display
 */
void display_cleanup(void) {
    if (sdlrenderer) SDL_DestroyRenderer(sdlrenderer);
    if (sdlwindow) SDL_DestroyWindow(sdlwindow);
    SDL_Quit();
}