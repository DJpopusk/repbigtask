#include "lodepng.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define radius 1
#define epsilon 50

typedef struct Node {
    unsigned char r, g, b, a;
    struct Node* up, * down, * left, * right, * parent;
    int rank;
} Node;

void gaussian_filter(unsigned char* image, unsigned w, unsigned h) {
    int x, y, i, j;
    float r, g, b;
    float kernel[3][3] = {
        {1.0 / 16, 1.0 / 8, 1.0 / 16},
        {1.0 / 8, 1.0 / 4, 1.0 / 8},
        {1.0 / 16, 1.0 / 8, 1.0 / 16}
    };
    unsigned char* tempImage = (unsigned char*)malloc(w * h * 4);

    for (y = 1; y < h - 1; y++) {
        for (x = 1; x < w - 1; x++) {
            r = g = b = 0;
            for (j = -radius; j <= radius; j++) {
                for (i = -radius; i <= radius; i++) {
                    int pixelPos = 4 * ((y + j) * w + (x + i));
                    r += image[pixelPos] * kernel[j + radius][i + radius];
                    g += image[pixelPos + 1] * kernel[j + radius][i + radius];
                    b += image[pixelPos + 2] * kernel[j + radius][i + radius];
                }
            }
            tempImage[4 * (y * w + x)] = (unsigned char)r;
            tempImage[4 * (y * w + x) + 1] = (unsigned char)g;
            tempImage[4 * (y * w + x) + 2] = (unsigned char)b;
            tempImage[4 * (y * w + x) + 3] = 255;
        }
    }

    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            image[4 * (y * w + x)] = tempImage[4 * (y * w + x)];
            image[4 * (y * w + x) + 1] = tempImage[4 * (y * w + x) + 1];
            image[4 * (y * w + x) + 2] = tempImage[4 * (y * w + x) + 2];
        }
    }

    free(tempImage);
}

void sobel_filter(unsigned char* image, unsigned w, unsigned h) {
    int x, y;
    int Gx[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}
    };
    int Gy[3][3] = {
        {1, 2, 1},
        {0, 0, 0},
        {-1, -2, -1}
    };
    unsigned char* tempImage = (unsigned char*)malloc(w * h * 4);

    for (y = 1; y < h - 1; y++) {
        for (x = 1; x < w - 1; x++) {
            float sumX = 0.0;
            float sumY = 0.0;
            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    unsigned char r = image[4 * ((y + i) * w + (x + j))];
                    sumX += Gx[i + 1][j + 1] * r;
                    sumY += Gy[i + 1][j + 1] * r;
                }
            }
            int sum = (int)sqrt(sumX * sumX + sumY * sumY);
            if (sum > 255) sum = 255;
            if (sum < 0) sum = 0;

            tempImage[4 * (y * w + x)] = sum;
            tempImage[4 * (y * w + x) + 1] = sum;
            tempImage[4 * (y * w + x) + 2] = sum;
            tempImage[4 * (y * w + x) + 3] = 255;
        }
    }

    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            image[4 * (y * w + x)] = tempImage[4 * (y * w + x)];
            image[4 * (y * w + x) + 1] = tempImage[4 * (y * w + x) + 1];
            image[4 * (y * w + x) + 2] = tempImage[4 * (y * w + x) + 2];
        }
    }

    free(tempImage);
}

void median_filter(unsigned char* image, unsigned w, unsigned h) {
    int windowSize = 3;
    int edge = windowSize / 2;
    unsigned char* tempImage = (unsigned char*)malloc(w * h * 4);

    for (int y = edge; y < h - edge; y++) {
        for (int x = edge; x < w - edge; x++) {
            unsigned char windowR[windowSize * windowSize];
            unsigned char windowG[windowSize * windowSize];
            unsigned char windowB[windowSize * windowSize];

            int k = 0;
            for (int i = -edge; i <= edge; i++) {
                for (int j = -edge; j <= edge; j++) {
                    int pixelPos = 4 * ((y + i) * w + (x + j));
                    windowR[k] = image[pixelPos];
                    windowG[k] = image[pixelPos + 1];
                    windowB[k] = image[pixelPos + 2];
                    k++;
                }
            }

            // Сортируем массивы
            for (int m = 0; m < windowSize * windowSize - 1; m++) {
                for (int n = m + 1; n < windowSize * windowSize; n++) {
                    if (windowR[m] > windowR[n]) {
                        unsigned char temp = windowR[m];
                        windowR[m] = windowR[n];
                        windowR[n] = temp;
                    }
                    if (windowG[m] > windowG[n]) {
                        unsigned char temp = windowG[m];
                        windowG[m] = windowG[n];
                        windowG[n] = temp;
                    }
                    if (windowB[m] > windowB[n]) {
                        unsigned char temp = windowB[m];
                        windowB[m] = windowB[n];
                        windowB[n] = temp;
                    }
                }
            }

            tempImage[4 * (y * w + x)] = windowR[(windowSize * windowSize) / 2];
            tempImage[4 * (y * w + x) + 1] = windowG[(windowSize * windowSize) / 2];
            tempImage[4 * (y * w + x) + 2] = windowB[(windowSize * windowSize) / 2];
            tempImage[4 * (y * w + x) + 3] = 255;
        }
    }

    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            image[4 * (y * w + x)] = tempImage[4 * (y * w + x)];
            image[4 * (y * w + x) + 1] = tempImage[4 * (y * w + x) + 1];
            image[4 * (y * w + x) + 2] = tempImage[4 * (y * w + x) + 2];
        }
    }

    free(tempImage);
}

int color_diff(Node* a, Node* b) {
    return sqrt(pow((a->r - b->r), 2) + pow((a->g - b->g), 2) + pow((a->b - b->b), 2));
}

Node* find(Node* node) {
    if (node->parent != node) {
        node->parent = find(node->parent);
    }
    return node->parent;
}

void union_sets(Node* a, Node* b) {
    Node* rootA = find(a);
    Node* rootB = find(b);

    if (rootA != rootB) {
        if (rootA->rank > rootB->rank) {
            rootB->parent = rootA;
        } else if (rootA->rank < rootB->rank) {
            rootA->parent = rootB;
        } else {
            rootB->parent = rootA;
            rootA->rank++;
        }
    }
}

Node* nodes;

void init_nodes(unsigned char* image, unsigned w, unsigned h) {
    nodes = (Node*)malloc(w * h * sizeof(Node));
    if (!nodes) {
        fprintf(stderr, "Malloc failed\n");
        exit(1);
    }

    for (unsigned y = 0; y < h; y++) {
        for (unsigned x = 0; x < w; x++) {
            Node* n = &nodes[y * w + x];
            n->r = image[4 * w * y + 4 * x + 0];
            n->g = image[4 * w * y + 4 * x + 1];
            n->b = image[4 * w * y + 4 * x + 2];
            n->a = image[4 * w * y + 4 * x + 3];
            n->up = (y == 0) ? NULL : &nodes[(y - 1) * w + x];
            n->down = (y == h - 1) ? NULL : &nodes[(y + 1) * w + x];
            n->left = (x == 0) ? NULL : &nodes[y * w + (x - 1)];
            n->right = (x == w - 1) ? NULL : &nodes[y * w + (x + 1)];
            n->parent = n;
            n->rank = 0;
        }
    }
}

void components(unsigned w, unsigned h) {
    for (unsigned y = 0; y < h; y++) {
        for (unsigned x = 0; x < w; x++) {
            Node* n = &nodes[y * w + x];

            // Skip white pixels (255, 255, 255)
            if (n->r == 0 && n->g == 0 && n->b == 0) continue;

            if (n->right && color_diff(n, n->right) < epsilon) {
                union_sets(n, n->right);
            }
            if (n->down && color_diff(n, n->down) < epsilon) {
                union_sets(n, n->down);
            }
        }
    }
}

void color_components(unsigned char* image, unsigned w, unsigned h) {
    srand(time(NULL));

    for (unsigned y = 0; y < h; y++) {
        for (unsigned x = 0; x < w; x++) {
            Node* n = &nodes[y * w + x];
            Node* root = find(n);

            // Skip black pixels (0, 0, 0)
            if (root->r <= 60 && root->g <= 60 && root->b <= 60) continue;

            // Assign random color to connected components
            if (root->parent == root) {
                root->r = rand() % 256;
                root->g = rand() % 256;
                root->b = rand() % 256;
            }

            image[4 * w * y + 4 * x + 0] = root->r;
            image[4 * w * y + 4 * x + 1] = root->g;
            image[4 * w * y + 4 * x + 2] = root->b;
            image[4 * w * y + 4 * x + 3] = 255; // Set alpha to opaque
        }
    }
}

void count_component_sizes(unsigned w, unsigned h, int* component_sizes) {
    for (unsigned y = 0; y < h; y++) {
        for (unsigned x = 0; x < w; x++) {
            Node* n = &nodes[y * w + x];
            Node* root = find(n);
            component_sizes[root - nodes]++;
        }
    }
}

void remove_small_components(unsigned char* image, unsigned w, unsigned h, int* component_sizes, int min_size) {
    for (unsigned y = 0; y < h; y++) {
        for (unsigned x = 0; x < w; x++) {
            Node* n = &nodes[y * w + x];
            Node* root = find(n);
            if (component_sizes[root - nodes] < min_size) {
                image[4 * w * y + 4 * x + 0] = 255; 
                image[4 * w * y + 4 * x + 1] = 255;
                image[4 * w * y + 4 * x + 2] = 255;
                image[4 * w * y + 4 * x + 3] = 255;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        printf("Please provide input PNG file\n");
        return 1;
    }

    unsigned error;
    unsigned char* image;
    unsigned width, height;
    error = lodepng_decode32_file(&image, &width, &height, argv[1]);
    if (error) {
        printf("Error %u: %s\n", error, lodepng_error_text(error));
        return 1;
    }

    gaussian_filter(image, width, height);
    printf("Применение фильтра Собеля\n");
    sobel_filter(image, width, height);
    printf("Применение медианного фильтра\n");
    median_filter(image, width, height);

    init_nodes(image, width, height);
    components(width, height);

    int* component_sizes = (int*)calloc(width * height, sizeof(int));
    count_component_sizes(width, height, component_sizes);

    remove_small_components(image, width, height, component_sizes, 50); // Минимальный размер компонента
    color_components(image, width, height);
    char outputFilename[] = "output.png";
    error = lodepng_encode32_file(outputFilename, image, width, height);
    if (error) {
        printf("Error %u: %s\n", error, lodepng_error_text(error));
    }

    free(image);
    free(nodes);
    free(component_sizes);
    return 0;
}