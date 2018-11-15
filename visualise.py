#!/usr/bin/env python3
import tkinter
import numpy as np
import random


class VisualScene:
    """
    Simple visual interface based on TKinter. Please replace with awesome visual interface.
    """
    color_list = ["yellow", "green", "cyan", "magenta", "red", "blue", "black", "brown"]
    directed_polygon = np.array([[0, -1], [1, 1], [-1, 1]])

    def __init__(self, data_reader):
        """
        Initializes a visual interface for the simulation. Updates every fixed amount of seconds.
        Represents the scene on a canvas.
        Important: this class progresses the simulation. After each drawing and potential delay,
        the visualisation calls for the progression to the next time step.
        Although it might be cleaner to move that to the simulation manager.
        :param data_reader: Scene to be drawn. The size of the scene is independent of the size of the visualization
        :return: None
        """
        self.data_reader = data_reader
        self.window = None
        self.env = None
        self.step_callback = None  # set in manager
        self.original_env = None
        self.canvas = None
        self.auto_loop = True
        self.p_size = 0.03
        self.colors = [random.choice(self.color_list) for _ in range(self.data_reader.num_particles)]
        self.time_delay = 100
        self.window = tkinter.Tk()
        self.window.title("Terrier")
        self._size = np.array((1000, 500))
        self.window.geometry("%dx%d" % (self._size[0], self._size[1]))
        self.window.grid()
        self.slice_num = 0

        self.canvas = tkinter.Canvas(self.window, bd=0, highlightthickness=0)
        self.canvas.grid(row=0, sticky=tkinter.W + tkinter.E + tkinter.N + tkinter.S)
        self.canvas.pack(fill=tkinter.BOTH, expand=1)
        self.window.bind("<Button-2>", self._give_position)

        # Geometry data
        x_centers = [-(self.data_reader.circle_radius + self.data_reader.circle_distance / 2),
                     (self.data_reader.circle_radius + self.data_reader.circle_distance / 2)]
        radii = [self.data_reader.circle_radius, self.data_reader.circle_radius]
        if self.data_reader.gate_radius > 0:
            x_centers.append(0)
            radii.append(self.data_reader.gate_radius)
        self.geom = tuple((x_centers, radii))
        if not self.auto_loop:
            self.disable_loop()
            print("Auto-updating of backend disabled. Press <Space> or click to advance simulation")
        # self.window.pack(fill=tkinter.BOTH,expand=1)

    @property
    def size(self):
        return [self.canvas.winfo_width(), self.canvas.winfo_height()]

    @size.setter
    def size(self, value):
        self._size = value
        self.window.geometry("%dx%d" % tuple(value))

    def start(self):
        """
        Starts the visualization loop
        :return:
        """
        self.loop()
        self.window.mainloop()

    def disable_loop(self):
        """
        Stop automatically redrawing.
        Enables space and Left-mouse click as progressing simulation.
        :return: None
        """
        self.auto_loop = False
        self.window.bind("<Button-1>", self.loop)
        self.window.bind("<space>", self.loop)

    def loop(self, _=None):
        """
        Public interface for visual scene loop. Has a callback reference to itself to keep the simulation going.
        :param _: Event object from tkinter
        :return: None
        """
        self.draw_scene()
        self.slice_num += 1
        self.data_reader.read_line()
        if self.auto_loop:
            self.window.after(self.time_delay, self.loop)

    def finish(self):
        """
        Cleanup. Called after self.start() returns

        :return:
        """
        self.window.destroy()
        self.auto_loop = False

    def _give_position(self, event):
        """
        Return the position of the click, reverted to the used coordinate system (in first quadrant of unit square)
        :param event: Mouse click or whatever
        :return: Coordinates relative to scene.
        """
        rel_x, rel_y = (event.x / self.size[0], 1 - event.y / self.size[1])
        x = (rel_x * 2 - 1) * self.data_reader.size[0]
        y = (rel_y * 2 - 1) * self.data_reader.size[1]
        print("Mouse location: (%.2f,%.2f)" % (x, y))

    def draw_scene(self):
        """
        Method that orders the draw commands of all objects within the scene.
        All objects are removed prior to the drawing step.
        :return: None
        """
        self.canvas.delete('all')
        self.draw_outline()
        self.draw_pedestrians()
        self.store_scene(None, 'images/slice-%04d.eps' % self.slice_num)

    def store_scene(self, _, filename=None):
        """
        Store a snapshot of the scene as an vector image.
        :param _: Event argument supplied by tkinter, can be ignored.
        :param filename: image file name. Leave empty for time-based (unique) name.

        :return: None
        """
        directory = 'images'
        if not filename:
            import time

            name = "scene#%d" % time.time()
            filename = "%s/%s-%.2f.eps" % (directory, name, self.data_reader.time)
        print("Snapshot at %.2f. Storing in %s" % (self.data_reader.time, filename))
        self.canvas.postscript(file=filename, pageheight=self.size[1], pagewidth=self.size[0])

    def draw_outline(self):
        for i in range(len(self.geom[0])):
            x_center = self.geom[0][i]
            radius = self.geom[1][i]
            rel_pos_array = np.array([x_center, 0]) / self.data_reader.size
            rel_size_array = radius / self.data_reader.size * self.size
            vis_pos_array = np.array([rel_pos_array[0] / 2 + 0.5, 1 - (rel_pos_array[1] / 2 + 0.5)]) * self.size
            start_pos_array = vis_pos_array - 0.5 * rel_size_array
            end_pos_array = vis_pos_array + 0.5 * rel_size_array
            self.canvas.create_oval(start_pos_array[0], start_pos_array[1], end_pos_array[0], end_pos_array[1], width=3)
        self.draw_bridge()

    def draw_pedestrians(self):
        """
        Draws all the pedestrians in the scene using the visual_pedestrian coordinates.
        :return: None
        """
        start_pos_array, end_pos_array = self.get_visual_pedestrian_coordinates()
        for index in range(self.data_reader.num_particles):
            self.canvas.create_oval(start_pos_array[index, 0], start_pos_array[index, 1],
                                    end_pos_array[index, 0], end_pos_array[index, 1], fill=self.colors[index])

    def get_visual_pedestrian_coordinates(self):
        """
        Computes the coordinates of all pedestrian relative to the visualization.
        Uses vectorized operations for speed increments
        :return: relative start coordinates, relative end coordinates.
        """
        rel_pos_array = self.data_reader.positions / self.data_reader.size
        rel_size_array = np.ones(
            self.data_reader.positions.shape) * self.p_size / self.data_reader.size * self.size
        vis_pos_array = np.hstack(
            (rel_pos_array[:, 0][:, None] / 2 + 0.5, 1 - (rel_pos_array[:, 1][:, None] / 2 + 0.5))) * self.size
        start_pos_array = vis_pos_array - 0.5 * rel_size_array
        end_pos_array = vis_pos_array + 0.5 * rel_size_array
        return start_pos_array, end_pos_array

    def draw_bridge(self):
        """
        Draws a line segment in the scene on its relative location. Can be used for (debugging) paths
        :return: None
        """
        x_0 = self.convert_relative_coordinate(
            np.array((-self.data_reader.bridge_size / 2, -self.data_reader.bridge_height / 2)) / self.data_reader.size)
        x_1 = self.convert_relative_coordinate(
            np.array((self.data_reader.bridge_size / 2, -self.data_reader.bridge_height / 2)) / self.data_reader.size)
        x_2 = self.convert_relative_coordinate(
            np.array((-self.data_reader.bridge_size / 2, self.data_reader.bridge_height / 2)) / self.data_reader.size)
        x_3 = self.convert_relative_coordinate(
            np.array((self.data_reader.bridge_size / 2, self.data_reader.bridge_height / 2)) / self.data_reader.size)
        self.canvas.create_rectangle(tuple(x_0) + tuple(x_3), fill='white', outline='white', width=3)
        self.canvas.create_line(tuple(x_0) + tuple(x_1), fill='black', width=2)
        self.canvas.create_line(tuple(x_2) + tuple(x_3), fill='black', width=2)

    def convert_relative_coordinate(self, coord):
        """
        Converts relative coordinates (from [-1,1]x[-1,1]) to screen size coordinates.
        Should raise an error when coordinates fall from scene,
        but this method is used so frequently I'd rather not make the computation
        Also changes the orientation to a Carthesian coordinate system
        :param coord: coordinates (fractions to be converted)
        :return: a Size with the coordinates of screen
        """
        return np.array([coord[0] / 2 + 0.5, 1 - (coord[1] / 2 + 0.5)]) * self.size


class DataReader:

    def __init__(self, filename=None):
        self.circle_radius = 1
        self.filename = 'results.dat'
        if filename is not None:
            self.filename = filename
        self.file = open(self.filename, 'r')
        self.num_particles = self.gate_radius = self.circle_radius = 0
        self.circle_distance = self.bridge_height = self.bridge_size = 0
        self.read_parameters()
        self.size = np.array([self.circle_radius * 4 + self.circle_distance, self.circle_radius * 2]) * 0.55
        self.times = []
        self.positions = np.zeros((self.num_particles, 2))
        self.directions = np.zeros(self.num_particles)
        self.counter = 0
        self.dt = 0.01
        self.time = 0

    def read_parameters(self):
        names = self.file.readline().split()
        parameters = self.file.readline().split()
        for name, parameter in zip(names, parameters):
            if not hasattr(self, name):
                print("Jammer, %s zit er niet in ", name)
            else:
                setattr(self, name, float(parameter))
        self.num_particles = int(self.num_particles)

    def read_line(self):
        self.time = float(self.file.readline().strip())
        self.times.append(self.time)
        raw_pos_x = self.file.readline().strip().split(' ')
        raw_pos_y = self.file.readline().strip().split(' ')
        raw_dirs = self.file.readline().strip().split(' ')
        for i in range(self.num_particles):
            self.positions[i, 0] = float(raw_pos_x[i])
            self.positions[i, 1] = float(raw_pos_y[i])
            self.directions[i] = float(raw_dirs[i])


if __name__ == '__main__':
    dr = DataReader()
    vs = VisualScene(dr)
    vs.start()
