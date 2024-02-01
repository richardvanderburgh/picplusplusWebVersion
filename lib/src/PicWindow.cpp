#include "PicWindow.h"
#include <iostream>

namespace UI {
    PicWindow::PicWindow() {
        if (SDL_Init(SDL_INIT_VIDEO) < 0)
        {
            std::cerr << "Failed to initialize the SDL2 library\n";
            return;
        }
        mp_window = std::unique_ptr <SDL_Window, SDLWindowDestroyer>(SDL_CreateWindow(WINDOW_NAME.c_str(),
            SDL_WINDOWPOS_CENTERED,
            SDL_WINDOWPOS_CENTERED,
            680, 480,
            0), SDLWindowDestroyer());
        if (!mp_window)
        {
            std::cerr << "Failed to create window\n";
            return;
        }
    }
	int PicWindow::render() {

        mp_surface = std::unique_ptr<SDL_Surface, SDLSurfaceDestroyer> (SDL_GetWindowSurface(mp_window.get()), SDLSurfaceDestroyer());

        if (!mp_surface)
        {
            std::cout << "Failed to get the surface from the window\n";
            return -1;
        }

        SDL_UpdateWindowSurface(mp_window.get());

        SDL_Delay(5000);
        return 1;
	}
}