# %%
import sys
import numpy as np
from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import *
import time
from PyQt5.QtGui import QPixmap, QIcon , QImage
from PyQt5 import QtCore, QtGui, QtWidgets
import time
from PyQt5.QtWidgets import QApplication, QSplashScreen
import ui
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from matplotlib.patches import Arc
import math
from io import BytesIO
import random
# Create a new figure
fig, ax = plt.subplots()

class ENP_ME(QMainWindow, ui.Ui_MainWindow):
    def __init__(self, motor_lenght, phases, conducters, courant, inner_rotor_radius, outer_rotor_radius, num_magnets, magnet_radius, num_ancoches, ancoche_radius, outer_ancoche_radius, outer_stator_radius, meshing_degrees, meshing_radius):
        super().__init__()
        self.setupUi(self)
        self.setWindowFlags(self.windowFlags())
        icon = QIcon("img\ENP Logo.png")
        self.setWindowIcon(icon)
        # pixmap = QPixmap("img/machine.png")
        # self.picture_place.setPixmap(pixmap)

        # Initialize UI fields with provided values
        self.inner_rotor_radius_ui.setText(str(inner_rotor_radius))
        self.outer_rotor_radius_ui.setText(str(outer_rotor_radius))
        self.num_magnets_ui.setText(str(num_magnets))
        self.magnet_radius_ui.setText(str(magnet_radius))
        self.num_ancoches_ui.setText(str(num_ancoches))
        self.ancoche_radius_ui.setText(str(ancoche_radius))
        self.outer_ancoche_radius_ui.setText(str(outer_ancoche_radius))
        self.outer_stator_radius_ui.setText(str(outer_stator_radius))
        self.meshing_degrees_ui.setText(str(meshing_degrees))
        self.meshing_radius_ui.setText(str(meshing_radius))
        self.motor_lenght_ui.setText(str(motor_lenght))
        self.phases_ui.setText(str(phases))
        self.conducters_ui.setText(str(conducters))
        self.courant_ui.setText(str(courant))

        self.update_parameters()

        self.draw_motor()
        # self.rayon()
        # self.points_permiabilite()
        # self.points_relectances()
        # self.pointsFMM()
        # self.potentiel()
        # self.flux()
        # self.heatmap()

        self.draw.clicked.connect(self.drawclocked)
        self.mesh.clicked.connect(self.meshclocked)
        self.calculate.clicked.connect(self.calculateclocked)

    def update_parameters(self):
        self.inner_rotor_radius = float(self.inner_rotor_radius_ui.text())
        self.outer_rotor_radius = float(self.outer_rotor_radius_ui.text())
        self.num_magnets = int(self.num_magnets_ui.text())
        self.magnet_radius = float(self.magnet_radius_ui.text())
        self.num_ancoches = int(self.num_ancoches_ui.text())
        self.ancoche_radius = float(self.ancoche_radius_ui.text())
        self.outer_ancoche_radius = float(self.outer_ancoche_radius_ui.text())
        self.outer_stator_radius = float(self.outer_stator_radius_ui.text())
        self.meshing_degrees = float(self.meshing_degrees_ui.text())
        self.meshing_radius = float(self.meshing_radius_ui.text())
        self.motor_lenght = float(self.motor_lenght_ui.text())
        self.phases = int(self.phases_ui.text())
        self.conducters = int(self.conducters_ui.text())
        self.courant = float(self.courant_ui.text())
        if self.meshing_radius > 0.001:
            self.meshing_radius = 0.001
            self.meshing_radius_ui.setText(str(self.meshing_radius))

    def drawclocked(self):
        self.picture_place.clear()
        self.update_parameters()
        self.draw_motor()
    def meshclocked(self):
        self.picture_place.clear()
        self.update_parameters()
        self.draw_motor()
        self.meshing()

    def calculateclocked(self):
        self.picture_place.clear()
        option = self.calculate_options.currentText()
        if option == "Potentiel":
            self.meshing()
            self.rayon()
            self.points_permiabilite()
            self.points_relectances()
            self.pointsFMM()
            self.potentiel()
            self.potentiel_heatmap()
        if option == "FMM":
            self.meshing()
            self.rayon()
            self.points_permiabilite()
            self.points_relectances()
            self.pointsFMM()
            self.FMM_heatmap()
        if option == "Flux":
            self.meshing()
            self.rayon()
            self.points_permiabilite()
            self.points_relectances()
            self.pointsFMM()
            self.potentiel()
            self.flux()
            self.flux_heatmap()
        if option == "Induction":
            self.meshing()
            self.rayon()
            self.points_permiabilite()
            self.points_relectances()
            self.pointsFMM()
            self.potentiel()
            self.flux()
            self.induction_heatmap()

    def draw_motor(self):
        fig, ax = plt.subplots()

        # Draw the inner rotor
        inner_rotor = patches.Circle((0, 0), radius=self.inner_rotor_radius, fill=False, color='b')
        ax.add_patch(inner_rotor)
        
        # Draw the outer rotor
        outer_rotor = patches.Circle((0, 0), radius=self.outer_rotor_radius, fill=False, color='b')
        ax.add_patch(outer_rotor)
        
        # Draw magnets
        angles = np.linspace(0, 2 * np.pi, 2 * self.num_magnets, endpoint=False)
        x = self.magnet_radius * np.cos(angles)
        y = self.magnet_radius * np.sin(angles)
        x1 = self.outer_rotor_radius * np.cos(angles)
        y1 = self.outer_rotor_radius * np.sin(angles)
        self.magnets_angles = angles

        for i in range(0, 2 * self.num_magnets - 1, 2):
            start_angle = np.degrees(np.arctan2(y[i], x[i]))
            end_angle = np.degrees(np.arctan2(y[i + 1], x[i + 1]))
            arc = Arc((0, 0), 2 * self.magnet_radius, 2 * self.magnet_radius, 
                      theta1=start_angle, theta2=end_angle, edgecolor='blue', linewidth=1)
            ax.add_patch(arc)
        
        for i in range(2 * self.num_magnets):
            ax.plot([x1[i], x[i]], [y1[i], y[i]], 'b')

        # Draw ancoches
        angles = np.linspace(0, 2 * np.pi, 2 * self.num_ancoches, endpoint=False)
        x2 = self.ancoche_radius * np.cos(angles)
        y2 = self.ancoche_radius * np.sin(angles)
        x3 = self.outer_ancoche_radius * np.cos(angles)
        y3 = self.outer_ancoche_radius * np.sin(angles)
        self.ancoche_angles = angles
        for i in range(0, 2 * self.num_ancoches - 1, 2):
            start_angle = np.degrees(np.arctan2(y2[i], x2[i]))
            end_angle = np.degrees(np.arctan2(y2[i + 1], x2[i + 1]))
            arc = Arc((0, 0), 2 * self.ancoche_radius, 2 * self.ancoche_radius, 
                      theta1=start_angle, theta2=end_angle, edgecolor='blue', linewidth=1)
            ax.add_patch(arc)
            arc = Arc((0, 0), 2 * self.outer_ancoche_radius, 2 * self.outer_ancoche_radius, 
                      theta1=start_angle + 180 / self.num_ancoches, theta2=end_angle + 180 / self.num_ancoches, 
                      edgecolor='blue', linewidth=1)
            ax.add_patch(arc)

        for i in range(2 * self.num_ancoches):
            ax.plot([x2[i], x3[i]], [y2[i], y3[i]], 'b')

        # Draw the outer stator
        outer_stator = patches.Circle((0, 0), radius=self.outer_stator_radius, fill=False, color='b')
        ax.add_patch(outer_stator)

        # Set the aspect of the plot to be equal
        ax.set_aspect('equal')
        ax.set_xlim(-0.070, 0.070)
        ax.set_ylim(-0.070, 0.070)
        ax.set_xlabel('X-axis')
        ax.set_ylabel('Y-axis')
        ax.axis('off')

        # Save the figure to a BytesIO object
        buf = BytesIO()
        fig.savefig(buf, format='png', bbox_inches='tight', pad_inches=0)
        buf.seek(0)

        # Convert the BytesIO object to a QImage
        qimage = QImage.fromData(buf.read())

        # Convert the QImage to a QPixmap
        pixmap = QPixmap.fromImage(qimage)

        # Set the QPixmap to the QLabel
        self.picture_place.setPixmap(pixmap)

        # # # Close the figure to free memory
        plt.close(fig)




    def meshing(self):
        
        self.fig, self.ax = plt.subplots()

        # Draw the inner rotor
        inner_rotor = patches.Circle((0, 0), radius=self.inner_rotor_radius, fill=False, color='b')
        self.ax.add_patch(inner_rotor)
        
        # Draw the outer rotor
        outer_rotor = patches.Circle((0, 0), radius=self.outer_rotor_radius, fill=False, color='b')
        self.ax.add_patch(outer_rotor)
        
        # Draw magnets
        angles = np.linspace(0, 2 * np.pi, 2 * self.num_magnets, endpoint=False)
        x = self.magnet_radius * np.cos(angles)
        y = self.magnet_radius * np.sin(angles)
        x1 = self.outer_rotor_radius * np.cos(angles)
        y1 = self.outer_rotor_radius * np.sin(angles)

        for i in range(0, 2 * self.num_magnets - 1, 2):
            start_angle = np.degrees(np.arctan2(y[i], x[i]))
            end_angle = np.degrees(np.arctan2(y[i + 1], x[i + 1]))
            arc = Arc((0, 0), 2 * self.magnet_radius, 2 * self.magnet_radius, 
                      theta1=start_angle, theta2=end_angle, edgecolor='blue', linewidth=1)
            self.ax.add_patch(arc)
        
        for i in range(2 * self.num_magnets):
            self.ax.plot([x1[i], x[i]], [y1[i], y[i]], 'b')

        # Draw ancoches
        angles = np.linspace(0, 2 * np.pi, 2 * self.num_ancoches, endpoint=False)
        x2 = self.ancoche_radius * np.cos(angles)
        y2 = self.ancoche_radius * np.sin(angles)
        x3 = self.outer_ancoche_radius * np.cos(angles)
        y3 = self.outer_ancoche_radius * np.sin(angles)

        for i in range(0, 2 * self.num_ancoches - 1, 2):
            start_angle = np.degrees(np.arctan2(y2[i], x2[i]))
            end_angle = np.degrees(np.arctan2(y2[i + 1], x2[i + 1]))
            arc = Arc((0, 0), 2 * self.ancoche_radius, 2 * self.ancoche_radius, 
                      theta1=start_angle, theta2=end_angle, edgecolor='blue', linewidth=1)
            self.ax.add_patch(arc)
            arc = Arc((0, 0), 2 * self.outer_ancoche_radius, 2 * self.outer_ancoche_radius, 
                      theta1=start_angle + 180 / self.num_ancoches, theta2=end_angle + 180 / self.num_ancoches, 
                      edgecolor='blue', linewidth=1)
            self.ax.add_patch(arc)

        for i in range(2 * self.num_ancoches):
            self.ax.plot([x2[i], x3[i]], [y2[i], y3[i]], 'b')

        # Draw the outer stator
        outer_stator = patches.Circle((0, 0), radius=self.outer_stator_radius, fill=False, color='b')
        self.ax.add_patch(outer_stator)

        # Parameters for the circle
        radius = float(self.meshing_radius)
        num_points = int(360 / float(self.meshing_degrees)+1) # Number of points to generate on the circle
        meshing_linewidth = 0.1
        # Generate angles uniformly spaced between 0 and 2π
        theta = np.linspace(0, 2 * np.pi, num_points)
        self.points_angle = theta

        # Initialize empty arrays for stacking
        X = np.array([]).reshape(0, num_points)  # No rows initially, `num_points` columns
        Y = np.array([]).reshape(0, num_points)

        # Loop to generate points for different radius multiples
        for i in range(0, int(self.outer_stator_radius/radius)+1):  # Note the loop starts from 1 to 5
            x = i * radius * np.cos(theta)  # Scaled by 'i'
            y = i * radius * np.sin(theta)
            X = np.vstack((X, x))  # Stack vertically
            Y = np.vstack((Y, y))

            self.points_x = X
            self.points_y = Y
        # Plot the points
        for i in range(X.shape[0]):  # Plot each row of X and Y
            plt.plot(X[i], Y[i], 'o', color='black', markersize=meshing_linewidth)  # Plot with circle markers and label
        
        
        for j in range(num_points-1):  # Iterate over each angle
            points_x = [X[i][j] for i in range(len(X))]  # X-coordinates for this angle across circles
            points_y = [Y[i][j] for i in range(len(X))]  # Y-coordinates for this angle across circles
            plt.plot(points_x, points_y,color='black', alpha=0.5) #
        # Add arcs between points with the same radius
        for i in range(X.shape[0]):  # Iterate over each circle
            # Draw arcs from the first to the last point
            start_angle = 0  # Degrees
            end_angle = 360  # Full circle
            radius_arc = i * radius  # The current radius for the circle
            arc = patches.Arc((0, 0), 2 * radius_arc, 2 * radius_arc, angle=0, theta1=start_angle, theta2=end_angle, linestyle='-', linewidth=meshing_linewidth, color='black')
            self.ax.add_patch(arc)
        elements = self.points_x.shape[0] * self.points_x.shape[1]

        # Set the aspect of the plot to be equal
        self.ax.set_aspect('equal')
        self.ax.set_xlim(-0.070, 0.070)
        self.ax.set_ylim(-0.070, 0.070)
        self.ax.set_xlabel('X-axis')
        self.ax.set_ylabel('Y-axis')
        self.ax.axis('off')
        # Save the figure to a BytesIO object
        buf = BytesIO()
        self.fig.savefig(buf, format='png', bbox_inches='tight', pad_inches=0)
        buf.seek(0)
        # Convert the BytesIO object to a QImage
        qimage = QImage.fromData(buf.read())

        # Convert the QImage to a QPixmap
        pixmap = QPixmap.fromImage(qimage)

        # Set the QPixmap to the QLabel
        self.picture_place.setPixmap(pixmap)

        # # # Close the figure to free memory
        # plt.close(fig)


        
    def rayon(self):

        self.rayon_radius = np.zeros((self.points_x.shape[0], self.points_x.shape[1]))
        for i in range(self.points_x.shape[0]):
            for j in range(self.points_x.shape[1]):
                self.rayon_radius[i][j] = math.sqrt(self.points_x[i][j]**2 + self.points_y[i][j]**2)
        self.rayon_radius = self.rayon_radius[:, :-1]
        self.points_angle = self.points_angle[:-1]

    def points_permiabilite(self):
        magnet_angle = 2*np.pi / (2 * self.num_magnets)
        ancoche_angle = 2*np.pi / (2 * self.num_ancoches)
        self.meshing_points_permiabilite = np.zeros((self.rayon_radius.shape[0], self.rayon_radius.shape[1]))
        for i in range(self.rayon_radius.shape[0]):
            for j in range(self.rayon_radius.shape[1]):
                if self.rayon_radius[i][j] <= self.inner_rotor_radius:
                    self.meshing_points_permiabilite[i][j] = 0
                elif self.rayon_radius[i][j] > self.inner_rotor_radius and self.rayon_radius[i][j] <= self.outer_rotor_radius:
                    self.meshing_points_permiabilite[i][j] = 1000
                elif self.rayon_radius[i][j] > self.outer_rotor_radius and self.rayon_radius[i][j] <= self.magnet_radius:
                    for num in range(self.magnets_angles.shape[0]):
                        if num % 2 == 0:
                            if self.points_angle[j] >= self.magnets_angles[num] and self.points_angle[j] < self.magnets_angles[num] + magnet_angle:
                                self.meshing_points_permiabilite[i][j] = 1
                            else:
                                self.meshing_points_permiabilite[i][j] = 1
                elif self.rayon_radius[i][j] > self.magnet_radius and self.rayon_radius[i][j] <= self.ancoche_radius:
                    self.meshing_points_permiabilite[i][j] = 1
                elif self.rayon_radius[i][j] > self.ancoche_radius and self.rayon_radius[i][j] <= self.outer_ancoche_radius:
                    for num in range(self.ancoche_angles.shape[0]):
                        if num % 2 == 0:
                            if self.points_angle[j] >= self.ancoche_angles[num] and self.points_angle[j] < self.ancoche_angles[num] + ancoche_angle:
                                self.meshing_points_permiabilite[i][j] = 1000
                            else:
                                self.meshing_points_permiabilite[i][j] = 1
                else:
                    self.meshing_points_permiabilite[i][j] = 1000
        np.set_printoptions(threshold=np.inf)  # Set printing threshold to infinity
        self.meshing_points_permiabilite = self.meshing_points_permiabilite * 4 * np.pi * 0.0000001

    def points_relectances(self):
        for i in range(self.rayon_radius.shape[0]-1):
            for j in range(self.rayon_radius.shape[1]):
                if self.rayon_radius[i][j] > self.inner_rotor_radius:
                    constant_radial     = (math.log(self.rayon_radius[i + 1][j] / self.rayon_radius[i][j])) / (self.motor_lenght * self.points_angle[1])
                    constant_tangentiel = self.points_angle[1] / (self.motor_lenght * math.log(self.rayon_radius[i + 1][j] / self.rayon_radius[i][j]))
                else:
                    constant_radial = 0
        self.points_relectances_radial      = constant_radial / self.meshing_points_permiabilite
        self.points_relectances_tangentiel  = constant_tangentiel / self.meshing_points_permiabilite
    
    def pointsFMM(self):
        self.meshing_points_FMM = np.zeros((self.rayon_radius.shape[0], self.rayon_radius.shape[1]))
        magnet_angle = 2*np.pi / (2 * self.num_magnets)
        k_matrix = np.zeros((1, self.magnets_angles.shape[0]))

        k = 1
        t = 0
        
        for i in range(k_matrix.shape[1]):
            k_matrix[0][i] = k
            t += 1
            if t == 2:
                k = - k
                t = 0

        for i in range(self.rayon_radius.shape[0]):
            for j in range(self.rayon_radius.shape[1]):
                if self.rayon_radius[i][j] > self.outer_rotor_radius and self.rayon_radius[i][j] <= self.magnet_radius:
                    for num in range(self.magnets_angles.shape[0]):
                        if num % 2 == 0:
                            
                            if self.points_angle[j] >= self.magnets_angles[num] and self.points_angle[j] < self.magnets_angles[num] + magnet_angle:
                                self.meshing_points_FMM[i][j] =  k_matrix[0][num] * 0.4 * (self.rayon_radius[i + 1][j] - self.rayon_radius[i][j]) / self.meshing_points_permiabilite[i][j]
    
        q = self.num_ancoches /(2 * self.phases * self.num_magnets / 2)
        ancoche_width = 360 / (2 * self.num_ancoches)
        num_points = int(ancoche_width / (self.points_angle[1] * 360 / (2 * np.pi)))
        e = np.zeros((int(q * self.phases),int(2 * np.pi / self.points_angle[1])))
        for i in range(1):
            for j in range(e.shape[1]):
                if (self.points_angle[j] * 360 / (2 * np.pi)) >= 0 and (self.points_angle[j] * 360 / (2 * np.pi)) < ancoche_width:
                    e[i][j] = - self.conducters * self.courant / 2
                elif (self.points_angle[j] * 360 / (2 * np.pi)) >= ancoche_width and (self.points_angle[j] * 360 / (2 * np.pi)) < 2 * ancoche_width:
                    e[i][j] = (self.conducters * self.courant / ancoche_width) * (self.points_angle[j] * 360 / (2 * np.pi)) - 3 * (self.conducters * self.courant / 2)
                elif (self.points_angle[j] * 360 / (2 * np.pi)) >= 2 *ancoche_width and (self.points_angle[j] * 360 / (2 * np.pi)) < 19 * ancoche_width:
                    e[i][j] =  self.conducters * self.courant / 2
                elif (self.points_angle[j] * 360 / (2 * np.pi)) >= 19 *ancoche_width and (self.points_angle[j] * 360 / (2 * np.pi)) < 20 * ancoche_width:
                    e[i][j] = -(self.conducters * self.courant / ancoche_width) * (self.points_angle[j] * 360 / (2 * np.pi)) + 39 * (self.conducters * self.courant / 2)
                elif (self.points_angle[j] * 360 / (2 * np.pi)) >= 20 *ancoche_width and (self.points_angle[j] * 360 / (2 * np.pi)) < 37 * ancoche_width:
                    e[i][j] = - self.conducters * self.courant / 2
                elif (self.points_angle[j] * 360 / (2 * np.pi)) >= 37 *ancoche_width and (self.points_angle[j] * 360 / (2 * np.pi)) < 38 * ancoche_width:
                    e[i][j] = (self.conducters * self.courant / ancoche_width) * (self.points_angle[j] * 360 / (2 * np.pi)) - 75 * (self.conducters * self.courant / 2)
                elif (self.points_angle[j] * 360 / (2 * np.pi)) >= 38 *ancoche_width and (self.points_angle[j] * 360 / (2 * np.pi)) < 55 * ancoche_width:
                    e[i][j] =  self.conducters * self.courant / 2
                elif (self.points_angle[j] * 360 / (2 * np.pi)) >= 55 *ancoche_width and (self.points_angle[j] * 360 / (2 * np.pi)) < 56 * ancoche_width:
                    e[i][j] = -(self.conducters * self.courant / ancoche_width) * (self.points_angle[j] * 360 / (2 * np.pi)) + 111 * (self.conducters * self.courant / 2)
                elif (self.points_angle[j] * 360 / (2 * np.pi)) >= 56 *ancoche_width and (self.points_angle[j] * 360 / (2 * np.pi)) < 72 * ancoche_width:
                    e[i][j] = - self.conducters * self.courant / 2
        k  = 1
        tt = 0
        for i in range(1,e.shape[0]):
            for j in range(e.shape[1]):
                e[i][j] = k * e[i-1][int(j-2*num_points)]
        
        for i in range(int(q),int(2*q)):
            for j in range(e.shape[1]):
                e[i][j] = -e[i][j]



        e_total = np.zeros((1, e.shape[1]))
        for i in range(e.shape[0]):
            for j in range(e.shape[1]):
                e_total[0][j] += e[i][j]
    
        self.dent_FMM    = np.zeros((self.rayon_radius.shape[0], self.rayon_radius.shape[1]))
        zh          = self.outer_ancoche_radius - self.ancoche_radius
        for i in range(self.rayon_radius.shape[0]):
            if self.rayon_radius[i][0] >= self.ancoche_radius and self.rayon_radius[i][0] < self.outer_ancoche_radius:
                for j in range (self.rayon_radius.shape[1]):
                    x = self.points_angle[j] * 360 / (2 * np.pi) // 5
                    if (x % 2 ) == 0:
                        self.dent_FMM[i][j] = e_total[0][j] * self.meshing_radius / zh
        # self.meshing_points_FMM +=  self.dent_FMM
        # print(self.meshing_points_FMM)
        
    def potentiel(self):
        d = 0
        self.elimination =0
        for i in range(self.rayon_radius.shape[0]):
            for j in range (self.rayon_radius.shape[1]):
                if self.rayon_radius[i][0] > self.inner_rotor_radius:
                    d += 1
                else:
                    self.elimination += 1
        
        
        potentiel_constant = np.zeros((d ,d))

        FMM_constant = np.zeros((d,1))
        ii = 0          
        for i in range(self.rayon_radius.shape[0]):
            if self.rayon_radius[i][0] > self.inner_rotor_radius :
                for j in range (self.rayon_radius.shape[1]):

                    if j == 0:
                        if i == self.rayon_radius.shape[0] -1:
                            a = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1]) + 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][-1])  
                            b = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1])
                            c = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][-1])
                            d = 0
                            e = -2 / (self.points_relectances_radial[i][j] +  self.points_relectances_radial[i-1][j])
                        elif i ==(self.elimination / self.rayon_radius.shape[1]):
                            a = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1]) + 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][-1])  
                            b = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1])
                            c = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][-1])
                            d = -2 / (self.points_relectances_radial[i][j] +  self.points_relectances_radial[i+1][j])
                            e = 0

                        else:
                            a = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1]) + 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][-1])
                            b = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1])
                            c = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][-1])
                            d = -2 / (self.points_relectances_radial[i][j] +  self.points_relectances_radial[i+1][j])
                            e = -2 / (self.points_relectances_radial[i][j] +  self.points_relectances_radial[i-1][j])

                        a_place = i * self.rayon_radius.shape[1] + j 
                        b_place = i * self.rayon_radius.shape[1] + j + 1 
                        c_place = i * self.rayon_radius.shape[1] + j + self.rayon_radius.shape[1] -1 
                        d_place = (i + 1) * self.rayon_radius.shape[1] + j
                        e_place = (i - 1) * self.rayon_radius.shape[1] + j

                    elif j == self.rayon_radius.shape[1]-1:
                        if i == self.rayon_radius.shape[0] -1:
                            a = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][0]) + 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1]) 
                            b = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][0])
                            c = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1])
                            d = 0
                            e = -2 / (self.points_relectances_radial[i][j] +  self.points_relectances_radial[i-1][j])
                        elif i ==(self.elimination / self.rayon_radius.shape[1]):
                            a = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][0]) + 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1]) 
                            b = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][0])
                            c = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1])
                            d = -2 / (self.points_relectances_radial[i][j] +  self.points_relectances_radial[i+1][j])
                            e = 0
                        else:
                            a = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][0]) + 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1]) 
                            b = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][0])
                            c = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1])
                            d = -2 / (self.points_relectances_radial[i][j] +  self.points_relectances_radial[i+1][j])
                            e = -2 / (self.points_relectances_radial[i][j] +  self.points_relectances_radial[i-1][j])

                        a_place = i * self.rayon_radius.shape[1] + j 
                        b_place = i * self.rayon_radius.shape[1]
                        c_place = i * self.rayon_radius.shape[1] + j - 1
                        d_place = (i + 1) * self.rayon_radius.shape[1] + j 
                        e_place = (i - 1) * self.rayon_radius.shape[1] + j 
                    else:
                        if i == self.rayon_radius.shape[0] -1:
                            a = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1]) + 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1]) 
                            b = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1])
                            c = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1])
                            d = 0
                            e = -2 / (self.points_relectances_radial[i][j] +  self.points_relectances_radial[i-1][j])
                        elif i ==(self.elimination / self.rayon_radius.shape[1]):
                            a = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1]) + 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1]) 
                            b = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1])
                            c = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1])
                            d = -2 / (self.points_relectances_radial[i][j] +  self.points_relectances_radial[i+1][j])
                            e = 0
                        else:
                            a = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1]) + 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1])
                            b = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1])
                            c = -2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1])
                            d = -2 / (self.points_relectances_radial[i][j] +  self.points_relectances_radial[i+1][j])
                            e = -2 / (self.points_relectances_radial[i][j] +  self.points_relectances_radial[i-1][j])

                        a_place = i * self.rayon_radius.shape[1] + j 
                        b_place = i * self.rayon_radius.shape[1] + j + 1
                        c_place = i * self.rayon_radius.shape[1] + j - 1 
                        d_place = (i + 1) * self.rayon_radius.shape[1] + j 
                        e_place = (i - 1) * self.rayon_radius.shape[1] + j 
                    

                    a_place = a_place - self.elimination
                    b_place = b_place - self.elimination
                    c_place = c_place - self.elimination
                    d_place = d_place - self.elimination
                    e_place = e_place - self.elimination
            
                    if i == self.rayon_radius.shape[0] -1:
                        potentiel_constant[ii][a_place] = a
                        potentiel_constant[ii][b_place] = b
                        potentiel_constant[ii][c_place] = c
                        potentiel_constant[ii][e_place] = e
                    elif i ==(self.elimination / self.rayon_radius.shape[1]):
                        potentiel_constant[ii][a_place] = a
                        potentiel_constant[ii][b_place] = b
                        potentiel_constant[ii][c_place] = c
                        potentiel_constant[ii][d_place] = d
                    else:
                        potentiel_constant[ii][a_place] = a
                        potentiel_constant[ii][b_place] = b
                        potentiel_constant[ii][c_place] = c
                        potentiel_constant[ii][d_place] = d
                        potentiel_constant[ii][e_place] = e

                    if i == self.rayon_radius.shape[0] -1:
                        FMM_constant[ii, 0] = (
                            self.meshing_points_FMM[i][j] * (1 / (self.points_relectances_radial[i][j] + self.points_relectances_radial[i-1][j]) 
                            )
                            + self.meshing_points_FMM[i-1][j] / (self.points_relectances_radial[i][j] + self.points_relectances_radial[i-1][j]) 
                        )
                    elif i ==(self.elimination / self.rayon_radius.shape[1]):
                        FMM_constant[ii, 0] = (
                            self.meshing_points_FMM[i][j] * ( 
                            - 1 / (self.points_relectances_radial[i][j] + self.points_relectances_radial[i+1][j])) 
                            - self.meshing_points_FMM[i+1][j] / (self.points_relectances_radial[i][j] + self.points_relectances_radial[i+1][j])
                        )
                    else:
                        FMM_constant[ii, 0] = (
                            self.meshing_points_FMM[i][j] * (1 / (self.points_relectances_radial[i][j] + self.points_relectances_radial[i-1][j]) 
                            - 1 / (self.points_relectances_radial[i][j] + self.points_relectances_radial[i+1][j]))
                            + self.meshing_points_FMM[i-1][j] / (self.points_relectances_radial[i][j] + self.points_relectances_radial[i-1][j]) 
                            - self.meshing_points_FMM[i+1][j] / (self.points_relectances_radial[i][j] + self.points_relectances_radial[i+1][j])
                        )
                    
                    ii +=1

        self.V = np.linalg.inv(potentiel_constant) @ FMM_constant
        self.points_v = np.zeros((self.rayon_radius.shape[0], self.rayon_radius.shape[1]))
        ii = 0
        for i in range(self.rayon_radius.shape[0]):
            if self.rayon_radius[i][0] > self.inner_rotor_radius :
                for j in range(self.rayon_radius.shape[1]):
                    self.points_v[i][j] = self.V[ii]
                    ii +=1

        # print(self.V)
    def flux(self):
        self.meshing_points_flux_left   = np.zeros((self.points_x.shape[0], self.points_x.shape[1] -1))
        self.meshing_points_flux_right  = np.zeros((self.points_x.shape[0], self.points_x.shape[1] -1))
        self.meshing_points_flux_up     = np.zeros((self.points_x.shape[0], self.points_x.shape[1] -1))
        self.meshing_points_flux_down   = np.zeros((self.points_x.shape[0], self.points_x.shape[1] -1))
        self.points_induction           = np.zeros((self.points_x.shape[0], self.points_x.shape[1] -1))



        for i in range(self.rayon_radius.shape[0]):
            if self.rayon_radius[i][0] > self.inner_rotor_radius :
                for j in range (self.rayon_radius.shape[1]):
                    
                    if j == self.rayon_radius.shape[1]-1:
                        if i == self.rayon_radius.shape[0] - 1:
                            left    = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][0]) * (self.points_v[i][j] - self.points_v[i][0])
                            right   = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1]) * (self.points_v[i][j] - self.points_v[i][j-1])
                            up      = 0
                            down    = 2/(self.points_relectances_radial[i][j] + self.points_relectances_radial[i-1][j]) * (self.points_v[i][j] - self.points_v[i-1][j] - 0.5 * (self.meshing_points_FMM[i][j] + self.meshing_points_FMM[i-1][j]))
    
                        elif i ==(self.elimination / self.rayon_radius.shape[1]):
                            left    = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][0]) * (self.points_v[i][j] - self.points_v[i][0])
                            right   = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1]) * (self.points_v[i][j] - self.points_v[i][j-1])
                            up      = 2/(self.points_relectances_radial[i][j] + self.points_relectances_radial[i+1][j]) * (self.points_v[i][j] - self.points_v[i+1][j] + 0.5 * (self.meshing_points_FMM[i][j] + self.meshing_points_FMM[i+1][j]))
                            down    = 0
                        else:
                            left    = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][0]) * (self.points_v[i][j] - self.points_v[i][0])
                            right   = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1]) * (self.points_v[i][j] - self.points_v[i][j-1])
                            up      = 2/(self.points_relectances_radial[i][j] + self.points_relectances_radial[i+1][j]) * (self.points_v[i][j] - self.points_v[i+1][j] + 0.5 * (self.meshing_points_FMM[i][j] + self.meshing_points_FMM[i+1][j]))
                            down    = 2/(self.points_relectances_radial[i][j] + self.points_relectances_radial[i-1][j]) * (self.points_v[i][j] - self.points_v[i-1][j] - 0.5 * (self.meshing_points_FMM[i][j] + self.meshing_points_FMM[i-1][j]))
                    else:
                        if i == self.rayon_radius.shape[0] - 1:
                            left    = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1]) * (self.points_v[i][j] - self.points_v[i][j+1])
                            right   = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1]) * (self.points_v[i][j] - self.points_v[i][j-1])
                            up      = 0
                            down    = 2/(self.points_relectances_radial[i][j] + self.points_relectances_radial[i-1][j]) * (self.points_v[i][j] - self.points_v[i-1][j] - 0.5 * (self.meshing_points_FMM[i][j] + self.meshing_points_FMM[i-1][j]))

                        elif i ==(self.elimination / self.rayon_radius.shape[1]):
                            left    = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1]) * (self.points_v[i][j] - self.points_v[i][j+1])
                            right   = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1]) * (self.points_v[i][j] - self.points_v[i][j-1])
                            up      = 2/(self.points_relectances_radial[i][j] + self.points_relectances_radial[i+1][j]) * (self.points_v[i][j] - self.points_v[i+1][j] + 0.5 * (self.meshing_points_FMM[i][j] + self.meshing_points_FMM[i+1][j]))
                            down    = 0
                        else:
                            left    = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j+1]) * (self.points_v[i][j] - self.points_v[i][j+1])
                            right   = 2/(self.points_relectances_tangentiel[i][j] + self.points_relectances_tangentiel[i][j-1]) * (self.points_v[i][j] - self.points_v[i][j-1])
                            up      = 2/(self.points_relectances_radial[i][j] + self.points_relectances_radial[i+1][j]) * (self.points_v[i][j] - self.points_v[i+1][j] + 0.5 * (self.meshing_points_FMM[i][j] + self.meshing_points_FMM[i+1][j]))
                            down    = 2/(self.points_relectances_radial[i][j] + self.points_relectances_radial[i-1][j]) * (self.points_v[i][j] - self.points_v[i-1][j] - 0.5 * (self.meshing_points_FMM[i][j] + self.meshing_points_FMM[i-1][j]))

                    self.meshing_points_flux_left[i][j]     = left
                    self.meshing_points_flux_right[i][j]    = right
                    self.meshing_points_flux_up[i][j]       = up
                    self.meshing_points_flux_down[i][j]     = down
                    if i == self.rayon_radius.shape[0] - 1:
                        self.points_induction[i][j] = 0
                    else:
                        self.points_induction[i][j] = math.sqrt((self.meshing_points_flux_down[i][j]/(self.rayon_radius[i][j] * self.meshing_degrees * self.motor_lenght)) ** 2 + (self.meshing_points_flux_right[i][j]/((self.rayon_radius[i+1][j] -self.rayon_radius[i][j] ) * self.motor_lenght)) ** 2)
        k = self.courant * (0.8 / 2.8) / self.points_induction.max()        

        q = self.num_ancoches /(2 * self.phases * self.num_magnets / 2)
        ancoche_width = 360 / (2 * self.num_ancoches)
        num_points = int(ancoche_width / (self.points_angle[1] * 360 / (2 * np.pi)))
        max_2 = 0
        for i in range(self.points_induction.shape[0]):
            if self.rayon_radius[i][0] > self.ancoche_radius and self.rayon_radius[i][0] <= self.outer_ancoche_radius:
                for j in range(self.points_induction.shape[1]):
                    n = (self.points_angle[j] * 360 / (2 * np.pi)) // ancoche_width
                    if n%2 != 0 :
                        self.points_induction[i][j] = 0
                        self.meshing_points_flux_left[i][j] = 0
                    else:
                        self.points_induction[i][j] = self.points_induction[i][j] * 5
                        self.meshing_points_flux_left[i][j] = self.meshing_points_flux_left[i][j] *5
                    
            # elif self.rayon_radius[i][0] > self.outer_ancoche_radius:
            if self.rayon_radius[i][0] > self.inner_rotor_radius and self.rayon_radius[i][0] < self.outer_rotor_radius:
                for j in range(self.points_induction.shape[1]):
                    self.points_induction[i][j] = self.points_induction[i][j] * (k*8)
                    self.meshing_points_flux_left[i][j] = self.meshing_points_flux_left[i][j]* 5 * k
            for i in range(self.points_induction.shape[0]):
            #     if self.rayon_radius[i][0] > self.outer_ancoche_radius:
                for j in range(self.points_induction.shape[1]):
                    if self.points_induction[i][j] > self.courant * (0.8 / 2.8) :
                        self.points_induction[i][j] = self.courant * (0.8 / 2.8) 

            #


                    
    def potentiel_heatmap(self):


        # Initialize the value matrix
        # val = np.zeros((self.points_relectances_radial.shape[0], self.points_relectances_radial.shape[1]))
        val = self.points_v
        # Calculate values based on radial and tangential reflectances
        # for i in range(self.points_relectances_radial.shape[0]):
        #     for j in range(self.points_relectances_radial.shape[1]):
        #         val[i][j] = math.sqrt(self.points_relectances_radial[i][j] ** 2 + self.points_relectances_tangentiel[i][j] ** 2)
  
        # Convert polar coordinates to Cartesian coordinates
        rayon = np.array(self.rayon_radius)
        theta = np.array(self.points_angle)

        X, Y = rayon * np.cos(theta), rayon * np.sin(theta)

        # Extend X and Y to the shape (n+1, n+1) for pcolormesh, handling wrap-around
        n_rows, n_cols = rayon.shape
        X_ext = np.zeros((n_rows + 1, n_cols + 1))
        Y_ext = np.zeros((n_rows + 1, n_cols + 1))
        X_ext[:-1, :-1] = X
        Y_ext[:-1, :-1] = Y

        # Repeat the last row and column to fill X_ext and Y_ext
        X_ext[-1, :-1] = X[-1, :]
        X_ext[:-1, -1] = X[:, 0]  # Wrap around column
        X_ext[-1, -1] = X[-1, 0]  # Wrap around corner

        Y_ext[-1, :-1] = Y[-1, :]
        Y_ext[:-1, -1] = Y[:, 0]  # Wrap around column
        Y_ext[-1, -1] = Y[-1, 0]  # Wrap around corner

        # Extend val to handle the wrap-around for pcolormesh
        val_ext = np.zeros((n_rows + 1, n_cols + 1))
        val_ext[:-1, :-1] = val
        val_ext[-1, :-1] = val[-1, :]
        val_ext[:-1, -1] = val[:, 0]  # Wrap around column
        val_ext[-1, -1] = val[-1, 0]  # Wrap around corner

        # Create the heatmap using pcolormesh with the 'jet' colormap and 50% opacity
        heatmap = self.ax.pcolormesh(X_ext, Y_ext, val_ext, cmap='jet', alpha=0.5, shading='auto')
        self.fig.colorbar(heatmap, ax=self.ax, orientation='vertical', label='Value')

        # Show the plot
        plt.show()
             # Save the figure to a BytesIO object
        buf = BytesIO()
        self.fig.savefig(buf, format='png', bbox_inches='tight', pad_inches=0)
        buf.seek(0)

        # Convert the BytesIO object to a QImage
        qimage = QImage.fromData(buf.read())

        # Convert the QImage to a QPixmap
        pixmap = QPixmap.fromImage(qimage)

        # Set the QPixmap to the QLabel
        self.picture_place.setPixmap(pixmap)

        # # # Close the figure to free memory
        plt.close(fig)



    def FMM_heatmap(self):


        # Initialize the value matrix
        # val = np.zeros((self.points_relectances_radial.shape[0], self.points_relectances_radial.shape[1]))
        val = self.meshing_points_FMM
        # Calculate values based on radial and tangential reflectances
        # for i in range(self.points_relectances_radial.shape[0]):
        #     for j in range(self.points_relectances_radial.shape[1]):
        #         val[i][j] = math.sqrt(self.points_relectances_radial[i][j] ** 2 + self.points_relectances_tangentiel[i][j] ** 2)
  
        # Convert polar coordinates to Cartesian coordinates
        rayon = np.array(self.rayon_radius)
        theta = np.array(self.points_angle)

        X, Y = rayon * np.cos(theta), rayon * np.sin(theta)

        # Extend X and Y to the shape (n+1, n+1) for pcolormesh, handling wrap-around
        n_rows, n_cols = rayon.shape
        X_ext = np.zeros((n_rows + 1, n_cols + 1))
        Y_ext = np.zeros((n_rows + 1, n_cols + 1))
        X_ext[:-1, :-1] = X
        Y_ext[:-1, :-1] = Y

        # Repeat the last row and column to fill X_ext and Y_ext
        X_ext[-1, :-1] = X[-1, :]
        X_ext[:-1, -1] = X[:, 0]  # Wrap around column
        X_ext[-1, -1] = X[-1, 0]  # Wrap around corner

        Y_ext[-1, :-1] = Y[-1, :]
        Y_ext[:-1, -1] = Y[:, 0]  # Wrap around column
        Y_ext[-1, -1] = Y[-1, 0]  # Wrap around corner

        # Extend val to handle the wrap-around for pcolormesh
        val_ext = np.zeros((n_rows + 1, n_cols + 1))
        val_ext[:-1, :-1] = val
        val_ext[-1, :-1] = val[-1, :]
        val_ext[:-1, -1] = val[:, 0]  # Wrap around column
        val_ext[-1, -1] = val[-1, 0]  # Wrap around corner

        # Create the heatmap using pcolormesh with the 'jet' colormap and 50% opacity
        heatmap = self.ax.pcolormesh(X_ext, Y_ext, val_ext, cmap='jet', alpha=0.5, shading='auto')
        self.fig.colorbar(heatmap, ax=self.ax, orientation='vertical', label='Value')

        # Show the plot
        plt.show()
             # Save the figure to a BytesIO object
        buf = BytesIO()
        self.fig.savefig(buf, format='png', bbox_inches='tight', pad_inches=0)
        buf.seek(0)

        # Convert the BytesIO object to a QImage
        qimage = QImage.fromData(buf.read())

        # Convert the QImage to a QPixmap
        pixmap = QPixmap.fromImage(qimage)

        # Set the QPixmap to the QLabel
        self.picture_place.setPixmap(pixmap)

        # # # Close the figure to free memory
        plt.close(fig)
    
    
    def flux_heatmap(self):
        # Initialize the value matrix
        # val = np.zeros((self.points_relectances_radial.shape[0], self.points_relectances_radial.shape[1]))
        val = self.meshing_points_flux_left
  
        # Convert polar coordinates to Cartesian coordinates
        rayon = np.array(self.rayon_radius)
        theta = np.array(self.points_angle)

        X, Y = rayon * np.cos(theta), rayon * np.sin(theta)

        # Extend X and Y to the shape (n+1, n+1) for pcolormesh, handling wrap-around
        n_rows, n_cols = rayon.shape
        X_ext = np.zeros((n_rows + 1, n_cols + 1))
        Y_ext = np.zeros((n_rows + 1, n_cols + 1))
        X_ext[:-1, :-1] = X
        Y_ext[:-1, :-1] = Y

        # Repeat the last row and column to fill X_ext and Y_ext
        X_ext[-1, :-1] = X[-1, :]
        X_ext[:-1, -1] = X[:, 0]  # Wrap around column
        X_ext[-1, -1] = X[-1, 0]  # Wrap around corner

        Y_ext[-1, :-1] = Y[-1, :]
        Y_ext[:-1, -1] = Y[:, 0]  # Wrap around column
        Y_ext[-1, -1] = Y[-1, 0]  # Wrap around corner

        # Extend val to handle the wrap-around for pcolormesh
        val_ext = np.zeros((n_rows + 1, n_cols + 1))
        val_ext[:-1, :-1] = val
        val_ext[-1, :-1] = val[-1, :]
        val_ext[:-1, -1] = val[:, 0]  # Wrap around column
        val_ext[-1, -1] = val[-1, 0]  # Wrap around corner

        # Create the heatmap using pcolormesh with the 'jet' colormap and 50% opacity
        heatmap = self.ax.pcolormesh(X_ext, Y_ext, val_ext, cmap='jet', alpha=0.5, shading='auto')
        self.fig.colorbar(heatmap, ax=self.ax, orientation='vertical', label='Value')

        # Show the plot
        plt.show()
             # Save the figure to a BytesIO object
        buf = BytesIO()
        self.fig.savefig(buf, format='png', bbox_inches='tight', pad_inches=0)
        buf.seek(0)

        # Convert the BytesIO object to a QImage
        qimage = QImage.fromData(buf.read())

        # Convert the QImage to a QPixmap
        pixmap = QPixmap.fromImage(qimage)

        # Set the QPixmap to the QLabel
        self.picture_place.setPixmap(pixmap)

        # # # Close the figure to free memory
        plt.close(fig)
    def induction_heatmap(self):
        # Initialize the value matrix
        # val = np.zeros((self.points_relectances_radial.shape[0], self.points_relectances_radial.shape[1]))
        val = self.points_induction
        # Calculate values based on radial and tangential reflectances
        # for i in range(self.points_relectances_radial.shape[0]):
        #     for j in range(self.points_relectances_radial.shape[1]):
        #         val[i][j] = math.sqrt(self.points_relectances_radial[i][j] ** 2 + self.points_relectances_tangentiel[i][j] ** 2)
  
        # Convert polar coordinates to Cartesian coordinates
        rayon = np.array(self.rayon_radius)
        theta = np.array(self.points_angle)

        X, Y = rayon * np.cos(theta), rayon * np.sin(theta)

        # Extend X and Y to the shape (n+1, n+1) for pcolormesh, handling wrap-around
        n_rows, n_cols = rayon.shape
        X_ext = np.zeros((n_rows + 1, n_cols + 1))
        Y_ext = np.zeros((n_rows + 1, n_cols + 1))
        X_ext[:-1, :-1] = X
        Y_ext[:-1, :-1] = Y

        # Repeat the last row and column to fill X_ext and Y_ext
        X_ext[-1, :-1] = X[-1, :]
        X_ext[:-1, -1] = X[:, 0]  # Wrap around column
        X_ext[-1, -1] = X[-1, 0]  # Wrap around corner

        Y_ext[-1, :-1] = Y[-1, :]
        Y_ext[:-1, -1] = Y[:, 0]  # Wrap around column
        Y_ext[-1, -1] = Y[-1, 0]  # Wrap around corner

        # Extend val to handle the wrap-around for pcolormesh
        val_ext = np.zeros((n_rows + 1, n_cols + 1))
        val_ext[:-1, :-1] = val
        val_ext[-1, :-1] = val[-1, :]
        val_ext[:-1, -1] = val[:, 0]  # Wrap around column
        val_ext[-1, -1] = val[-1, 0]  # Wrap around corner

        # Create the heatmap using pcolormesh with the 'jet' colormap and 50% opacity
        heatmap = self.ax.pcolormesh(X_ext, Y_ext, val_ext, cmap='jet', alpha=0.5, shading='auto')
        self.fig.colorbar(heatmap, ax=self.ax, orientation='vertical', label='Value')

        # Show the plot
        plt.show()
             # Save the figure to a BytesIO object
        buf = BytesIO()
        self.fig.savefig(buf, format='png', bbox_inches='tight', pad_inches=0)
        buf.seek(0)

        # Convert the BytesIO object to a QImage
        qimage = QImage.fromData(buf.read())

        # Convert the QImage to a QPixmap
        pixmap = QPixmap.fromImage(qimage)

        # Set the QPixmap to the QLabel
        self.picture_place.setPixmap(pixmap)

        # # # Close the figure to free memory
        plt.close(fig)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = ENP_ME(
        motor_lenght=0.110, phases=3, conducters=34, courant=2.8, inner_rotor_radius=0.0198, 
        outer_rotor_radius=0.03625, num_magnets=4, magnet_radius=0.04025, num_ancoches=36, 
        ancoche_radius=0.04125, outer_ancoche_radius=0.05385, outer_stator_radius=0.06915, 
        meshing_degrees=3, meshing_radius=0.001
    )
    

img_splash = QSplashScreen(QPixmap("img\ENP Logo.png"))
img_splash.resize(400,350)
img_splash.show()
time.sleep(3)
img_splash.finish(window)

window.show()
sys.exit(app.exec_())


# app.exec_()



