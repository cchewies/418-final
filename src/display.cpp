#include "compact_defines.h"
#include <SDL2/SDL.h>
#include <iostream>

// Window size
int gWidth = 800;
int gHeight = 600;

// SDL objects
SDL_Window* gWindow = nullptr;
SDL_Renderer* gRenderer = nullptr;
SDL_Texture* gTexture = nullptr;
uint32_t* gBuffer = nullptr;

/**
 * @brief Handle keyboard input
 * 
 * @param event SDL2 event
 */
static void handle_key(const SDL_Event& event) {
    if (event.type == SDL_KEYDOWN && event.key.keysym.sym == SDLK_ESCAPE) {
        fprintf(stderr, "Quit");
        exit(0);
    }
}

/**
 * @brief Render a frame
 * 
 */
void display_render() {

    // Handle key events
    SDL_Event event;
    while (SDL_PollEvent(&event)) {
        handle_key(event);
    }

    // TEMP: Fill the buffer with a gradient
    for (int y = 0; y < gHeight; ++y) {
        for (int x = 0; x < gWidth; ++x) {
            gBuffer[y * gWidth + x] = (255 << 24) | (x * 255 / gWidth << 16) | (y * 255 / gHeight << 8);
        }
    }

    // Render stuff
    SDL_UpdateTexture(gTexture, nullptr, gBuffer, gWidth * sizeof(uint32_t));
    SDL_RenderClear(gRenderer);
    SDL_RenderCopy(gRenderer, gTexture, nullptr, nullptr);
    SDL_RenderPresent(gRenderer);
}

/**
 * @brief Initialize display
 */
void display_init(void) {
    SDL_Init(SDL_INIT_VIDEO);

    gWindow = SDL_CreateWindow("barnes hut",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
        gWidth, gHeight, 0);

    gRenderer = SDL_CreateRenderer(gWindow, -1, SDL_RENDERER_ACCELERATED);

    gTexture = SDL_CreateTexture(gRenderer,
        SDL_PIXELFORMAT_RGBA8888,
        SDL_TEXTUREACCESS_STREAMING,
        gWidth, gHeight);

    gBuffer = new uint32_t[gWidth * gHeight];
}
