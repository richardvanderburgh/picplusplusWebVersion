#include <SDL2/SDL.h>
#include <string>
#include <memory>
namespace UI {
    struct SDLWindowDestroyer
    {
        void operator()(SDL_Window* w) const
        {
            SDL_DestroyWindow(w);
        }
    };

    struct SDLSurfaceDestroyer
    {
        void operator()(SDL_Surface* w) const
        {
            SDL_FreeSurface(w);
        }
    };
    class Window {
    public:
        Window();
        int render();
    private:
        const std::string WINDOW_NAME = "PIC++";
        std::unique_ptr<SDL_Window, SDLWindowDestroyer> mp_window;
        //// Will need it's own class eventually
        std::unique_ptr<SDL_Surface, SDLSurfaceDestroyer> mp_surface;
    };
}