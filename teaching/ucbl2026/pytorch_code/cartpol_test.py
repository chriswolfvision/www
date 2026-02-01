import gymnasium as gym
import pygame
import time
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="pygame")

# Initialize environment
env = gym.make("CartPole-v1", render_mode="human")
env.reset()

print("Left/Right arrows; 'ESC' to quit.")

running = True
while running:
    action = None
    keys = pygame.key.get_pressed()

    if keys[pygame.K_LEFT]:
        action = 0
    elif keys[pygame.K_RIGHT]:
        action = 1

    if action is not None:
        obs, reward, term, succ, info = env.step(action)

        if term:
            print("\n" + "="*20 + "\n  EPISODE FAILED!\n" + "="*20 + "\n" )

            # Keep the window open for a second, then reset
            time.sleep(1.0)
            env.reset()

        elif succ:
            print("Time Limit Reached (Success!)")
            env.reset()

    for event in pygame.event.get():
        if event.type == pygame.QUIT or (event.type == pygame.KEYDOWN and event.key == pygame.K_ESCAPE):
            running = False

env.close()
pygame.quit()