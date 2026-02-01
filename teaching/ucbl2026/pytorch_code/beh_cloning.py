import gymnasium as gym
import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import random
from collections import deque
import time
import argparse

# Hyperparameters
DEVICE = torch.device("cuda" if torch.cuda.is_available() else "cpu")
LR = 1e-3
GAMMA = 0.99
MEM_SIZE = 10000
BATCH_SIZE = 64

class StudentNet(nn.Module):
    def __init__(self, state_dim, action_dim):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(state_dim, 32),
            nn.ReLU(),
            nn.Linear(32, action_dim)
        )
    def forward(self, x):
        return self.net(x)

def visualize_agent(model, episodes=2, max_steps_per_episode=30, env_name="CartPole-v1"):
    # Create a new environment instance with rendering enabled
    render_env = gym.make(env_name, render_mode="human")

    print(f"Visualizing {episodes} episodes for max len {max_steps_per_episode}...")

    for ep in range(episodes):
        state, _ = render_env.reset()
        done = False
        truncated = False
        score = 0

        stepnr=0
        while not (done or truncated):
            # Slow down the loop slightly so we can actually see the movement
            render_env.render()
            time.sleep(0.02)

            # Agent selects action
            state_t = torch.tensor(state, dtype=torch.float32).to(DEVICE)
            with torch.no_grad():
                action = torch.argmax(model(state_t)).item()

            state, reward, done, truncated, _ = render_env.step(action)
            score += reward
            stepnr += 1
            if stepnr>max_steps_per_episode:
                break

        print(f"Episode {ep + 1} finished with score: {score}")

    render_env.close()

def behavior_cloning(exp_states, exp_act, num_epochs):
    States=4
    Actions=2
    model = StudentNet(States, Actions).to(DEVICE)
    optimizer = optim.Adam(model.parameters(), lr=1e-3)
    criterion = nn.CrossEntropyLoss()

    # For simplicity, we do not batch:
    # do one forward pass with the full data
    for epoch in range(num_epochs):
        pred_act = model(exp_states)
        loss = criterion(pred_act, exp_act)

        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if epoch % 20 == 0:
            print(f"Epoch {epoch}, Loss: {loss.item():.4f}")
    print(f"Epoch {epoch}, Loss: {loss.item():.4f}")

    return model

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--rl",
        action="store_true",
        help="Train the RL agent, generate trajectories and store them",
    )
    parser.add_argument(
        "--load",
        action="store_true",
        help="Load trajectories and store them",
    )
    parser.add_argument(
        "--bc",
        action="store_true",
        help="Run behavior cloning",
    )
    parser.add_argument(
        "--max_len",
        type=int,
        default=100,
        help="Max Length of one episode when visualizing",
    )
    parser.add_argument(
        "--bc_epochs",
        type=int,
        default=100,
        help="Num epochs for behavior cloning",
    )
    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()

    if args.rl:
        print ("Code for teacher training not provided ...")

    if args.load:
        print ("Loading trajectories from file.")
        d=np.load("traj.npz")
        exp_states=d["exp_states"]
        exp_act=d["exp_act"]

    if args.bc:
        print ("Training student with behavior cloning.")
        exp_states = torch.tensor(exp_states, dtype=torch.float32).to(DEVICE)
        exp_act = torch.tensor(exp_act, dtype=torch.long).to(DEVICE)
        student = behavior_cloning(exp_states, exp_act, args.bc_epochs)
        print ("Student trained, visualizing...")
        visualize_agent(student, max_steps_per_episode=args.max_len)

    print ("Done.")