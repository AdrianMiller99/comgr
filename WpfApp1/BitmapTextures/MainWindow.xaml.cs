using System;
using System.Collections.Generic;
using System.Numerics;
using System.Windows;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace BitmapTextures;

/// <summary>
/// Loads and samples bitmap textures with gamma-correct conversion.
/// Converts sRGB textures to linear RGB for physically correct rendering.
/// </summary>
public class Texture
{
    private Vector3[] _pixels;
    private int _width;
    private int _height;

    /// <summary>
    /// Loads a texture from a file path.
    /// </summary>
    /// <param name="filePath">Path to the texture file (jpg, png, etc.).</param>
    public Texture(string filePath)
    {
        LoadFromFile(filePath);
    }

    /// <summary>
    /// Loads texture data from a file and converts from sRGB to linear RGB.
    /// </summary>
    /// <param name="filePath">Path to the texture file.</param>
    private void LoadFromFile(string filePath)
    {
        BitmapImage bitmapImage = new BitmapImage(new Uri(filePath, UriKind.RelativeOrAbsolute));
        
        // Convert to a format we can read (BGRA32)
        FormatConvertedBitmap convertedBitmap = new FormatConvertedBitmap();
        convertedBitmap.BeginInit();
        convertedBitmap.Source = bitmapImage;
        convertedBitmap.DestinationFormat = PixelFormats.Bgra32;
        convertedBitmap.EndInit();

        _width = convertedBitmap.PixelWidth;
        _height = convertedBitmap.PixelHeight;
        
        int stride = _width * 4;
        byte[] pixelData = new byte[_height * stride];
        convertedBitmap.CopyPixels(pixelData, stride, 0);

        // Convert from sRGB bytes to linear RGB float3
        _pixels = new Vector3[_width * _height];
        for (int y = 0; y < _height; y++)
        {
            for (int x = 0; x < _width; x++)
            {
                int idx = y * stride + x * 4;
                byte b = pixelData[idx + 0];
                byte g = pixelData[idx + 1];
                byte r = pixelData[idx + 2];
                
                // Convert sRGB to linear RGB (gamma correction)
                _pixels[y * _width + x] = new Vector3(
                    SrgbToLinear(r / 255f),
                    SrgbToLinear(g / 255f),
                    SrgbToLinear(b / 255f)
                );
            }
        }
    }

    /// <summary>
    /// Converts a single sRGB color component to linear RGB.
    /// </summary>
    /// <param name="srgb">sRGB value in range [0, 1].</param>
    /// <returns>Linear RGB value.</returns>
    /// <remarks>
    /// Uses the standard sRGB to linear conversion:
    /// - For dark values (≤ 0.04045): simple division
    /// - For others: power function with gamma 2.4
    /// </remarks>
    private static float SrgbToLinear(float srgb)
    {
        if (srgb <= 0.04045f)
            return srgb / 12.92f;
        else
            return MathF.Pow((srgb + 0.055f) / 1.055f, 2.4f);
    }

    /// <summary>
    /// Samples the texture at UV coordinates with bilinear filtering.
    /// </summary>
    /// <param name="u">Horizontal texture coordinate (wraps around).</param>
    /// <param name="v">Vertical texture coordinate (wraps around).</param>
    /// <returns>Linear RGB color sampled from the texture.</returns>
    public Vector3 Sample(float u, float v)
    {
        // Wrap UV coordinates to [0, 1] range
        u = u - MathF.Floor(u);
        v = v - MathF.Floor(v);
        
        // Convert to pixel coordinates
        float x = u * (_width - 1);
        float y = v * (_height - 1);
        
        // Bilinear interpolation
        int x0 = (int)MathF.Floor(x);
        int y0 = (int)MathF.Floor(y);
        int x1 = Math.Min(x0 + 1, _width - 1);
        int y1 = Math.Min(y0 + 1, _height - 1);
        
        float fx = x - x0;
        float fy = y - y0;
        
        Vector3 c00 = _pixels[y0 * _width + x0];
        Vector3 c10 = _pixels[y0 * _width + x1];
        Vector3 c01 = _pixels[y1 * _width + x0];
        Vector3 c11 = _pixels[y1 * _width + x1];
        
        // Bilinear blend
        Vector3 c0 = Vector3.Lerp(c00, c10, fx);
        Vector3 c1 = Vector3.Lerp(c01, c11, fx);
        
        return Vector3.Lerp(c0, c1, fy);
    }
}

/// <summary>
/// Represents a ray in 3D space for path tracing.
/// </summary>
public struct Ray
{
    public Vector3 E; // Origin
    public Vector3 d; // Direction
}

/// <summary>
/// Material properties for surfaces.
/// Supports solid colors or textures, with optional emission and roughness.
/// </summary>
public struct Material
{
    /// <summary>
    /// Base color (used if Texture is null, or multiplied with texture).
    /// </summary>
    public Vector3 Color;
    
    /// <summary>
    /// Light emission.
    /// </summary>
    public Vector3 Emission;
    
    /// <summary>
    /// Surface roughness: 0 = mirror, 1 = diffuse.
    /// </summary>
    public float Roughness;
    
    /// <summary>
    /// Optional texture map (null for solid color).
    /// </summary>
    public Texture Texture;
}

/// <summary>
/// Ray-surface intersection information.
/// </summary>
public struct HitInfo
{
    public bool DidHit;
    public float distance; // Distance along ray
    public Vector3 H; // Hit point
    public Vector3 n; // Normal
    public Vector2 UV; // Texture coordinates
    public Material Material;
}

/// <summary>
/// Sphere primitive with texture coordinate generation.
/// </summary>
public class Sphere
{
    public Vector3 C; // Center
    public float r; // Radius
    public Material Material;

    /// <summary>
    /// Computes ray-sphere intersection with UV coordinates.
    /// </summary>
    /// <param name="ray">The ray to test for intersection.</param>
    /// <returns>HitInfo with intersection details, or DidHit=false if no intersection.</returns>
    public HitInfo Intersect(Ray ray)
    {
        Vector3 OC = ray.E - C;
        float a = Vector3.Dot(ray.d, ray.d);
        float b = 2.0f * Vector3.Dot(OC, ray.d);
        float c = Vector3.Dot(OC, OC) - r * r;

        float discriminant = b * b - 4.0f * a * c;
        if (discriminant < 0.0f)
            return new HitInfo { DidHit = false };

        float sqrtDiscriminant = MathF.Sqrt(discriminant);
        float t0 = (-b - sqrtDiscriminant) / (2.0f * a);
        float t1 = (-b + sqrtDiscriminant) / (2.0f * a);

        const float EPS = 1e-3f;
        float tHit = t0;
        if (tHit < EPS)
        {
            tHit = t1;
            if (tHit < EPS)
                return new HitInfo { DidHit = false };
        }

        Vector3 H = ray.E + tHit * ray.d;
        Vector3 n = Vector3.Normalize(H - C);

        // Calculate UV coordinates from sphere normal
        Vector2 uv = CalculateSphereUV(n);

        return new HitInfo
        {
            DidHit = true,
            distance = tHit,
            H = H,
            n = n,
            UV = uv,
            Material = this.Material
        };
    }

    /// <summary>
    /// Calculates UV coordinates from a sphere's normal vector.
    /// </summary>
    /// <param name="normal">Normalized surface normal.</param>
    /// <returns>UV coordinates in range [0, 1].</returns>
    private Vector2 CalculateSphereUV(Vector3 normal)
    {
        // Spherical coordinates: theta (elevation), phi (azimuth)
        float theta = MathF.Acos(Math.Clamp(normal.Y, -1f, 1f));
        float phi = MathF.Atan2(normal.Z, normal.X);
        
        // Convert to UV (0 to 1 range)
        float u = 0.5f + phi / (2.0f * MathF.PI);
        float v = theta / MathF.PI;
        
        return new Vector2(u, v);
    }
}

/// <summary>
/// Scene container for all geometric objects.
/// </summary>
public class Scene
{
    public List<Sphere> Spheres = new List<Sphere>();

    /// <summary>
    /// Finds the closest intersection of a ray with any object in the scene.
    /// </summary>
    /// <param name="ray">The ray to trace through the scene.</param>
    /// <returns>HitInfo for the closest intersection, or DidHit=false if no intersection.</returns>
    public HitInfo Intersect(Ray ray)
    {
        HitInfo closest = new HitInfo { DidHit = false, distance = float.MaxValue };
        foreach (var s in Spheres)
        {
            var hit = s.Intersect(ray);
            if (hit.DidHit && hit.distance < closest.distance)
                closest = hit;
        }
        return closest;
    }
}

/// <summary>
/// Implements Monte Carlo path tracing with importance sampling and Russian Roulette.
/// </summary>
/// <remarks>
/// Path tracing simulates global illumination by recursively tracing light paths.
/// Key techniques used:
/// - Monte Carlo integration for rendering equation
/// - Cosine-weighted hemisphere sampling for diffuse surfaces
/// - Russian Roulette for unbiased path termination
/// </remarks>
public class PathTracer
{
    private Random _rng;
    private const int MinBounces = 3;
    private const float RRProbability = 0.8f;

    public PathTracer(int seed = 42)
    {
        _rng = new Random(seed);
    }

    /// <summary>
    /// Traces a ray through the scene and computes the incoming radiance.
    /// Uses Monte Carlo integration with Russian Roulette path termination.
    /// </summary>
    /// <param name="ray">The ray to trace.</param>
    /// <param name="scene">The scene to trace through.</param>
    /// <param name="depth">Current recursion depth.</param>
    /// <returns>The radiance (color) arriving along the ray.</returns>
    public Vector3 Trace(Ray ray, Scene scene, int depth = 0)
    {
        HitInfo hit = scene.Intersect(ray);
        if (!hit.DidHit)
            return Vector3.Zero; // Background is black

        Vector3 color = hit.Material.Emission;

        // Russian Roulette: probabilistically terminate paths after MinBounces
        if (depth >= MinBounces)
        {
            float rand = (float)_rng.NextDouble();
            if (rand > RRProbability)
                return color; // Terminate path here
        }

        // Compute reflected radiance using BRDF sampling
        Vector3 reflected = SampleBRDF(ray, hit, scene, depth);
        
        // Apply Russian Roulette compensation to keep estimator unbiased
        if (depth >= MinBounces)
            reflected /= RRProbability;
        
        color += reflected;

        return color;
    }

    /// <summary>
    /// Samples the BRDF and computes reflected radiance.
    /// </summary>
    /// <param name="ray">The incoming ray.</param>
    /// <param name="hit">The surface hit information.</param>
    /// <param name="scene">The scene.</param>
    /// <param name="depth">Current recursion depth.</param>
    /// <returns>The reflected radiance from this surface.</returns>
    private Vector3 SampleBRDF(Ray ray, HitInfo hit, Scene scene, int depth)
    {
        float roughness = hit.Material.Roughness;
        
        // Get surface color from texture or solid color
        Vector3 albedo = hit.Material.Texture != null
            ? hit.Material.Texture.Sample(hit.UV.X, hit.UV.Y)
            : hit.Material.Color;
        
        Vector3 d_r = Vector3.Normalize(ray.d - 2.0f * Vector3.Dot(ray.d, hit.n) * hit.n);
        
        Vector3 sampleDir;
        
        if (roughness < 0.5f)
        {
            // Reflective: importance sample around reflection
            if (roughness < 0.01f)
                sampleDir = d_r;
            else
                sampleDir = PerturbReflection(d_r, roughness);
        }
        else
        {
            // Diffuse: hemisphere sampling
            sampleDir = SampleHemisphereCosine(hit.n);
        }
        
        Ray scattered = new Ray 
        { 
            E = hit.H + hit.n * 1e-4f,
            d = sampleDir 
        };

        Vector3 incomingLight = Trace(scattered, scene, depth + 1);

        return albedo * incomingLight;
    }
    
    /// <summary>
    /// Perturbs a reflection direction for glossy materials.
    /// Creates a cone of directions around the perfect reflection.
    /// </summary>
    /// <param name="reflectionDir">Perfect reflection direction.</param>
    /// <param name="roughness">Surface roughness (controls cone width).</param>
    /// <returns>Perturbed reflection direction.</returns>
    private Vector3 PerturbReflection(Vector3 reflectionDir, float roughness)
    {
        Vector3 v = SampleHemisphereCosine(reflectionDir);
        float t = roughness * roughness;
        return Vector3.Normalize(Vector3.Lerp(reflectionDir, v, t));
    }
    
    /// <summary>
    /// Cosine-weighted hemisphere sampling for diffuse surfaces.
    /// </summary>
    private Vector3 SampleHemisphereCosine(Vector3 normal)
    {
        float r1 = (float)_rng.NextDouble();
        float r2 = (float)_rng.NextDouble();

        float phi = 2.0f * MathF.PI * r1;
        float sinTheta = MathF.Sqrt(r2);
        float cosTheta = MathF.Sqrt(1.0f - r2);

        Vector3 tangent = GetPerpendicular(normal);
        Vector3 bitangent = Vector3.Cross(normal, tangent);

        Vector3 direction = 
            tangent * (MathF.Cos(phi) * sinTheta) +
            bitangent * (MathF.Sin(phi) * sinTheta) +
            normal * cosTheta;

        return Vector3.Normalize(direction);
    }

    /// <summary>
    /// Gets an arbitrary perpendicular vector.
    /// </summary>
    private Vector3 GetPerpendicular(Vector3 n)
    {
        Vector3 a = MathF.Abs(n.X) > 0.9f 
            ? new Vector3(0, 1, 0) 
            : new Vector3(1, 0, 0);
        return Vector3.Normalize(Vector3.Cross(n, a));
    }
}

/// <summary>
/// Main window for the Bitmap Textures renderer.
/// Renders a path-traced scene with gamma-correct textured spheres.
/// </summary>
public partial class MainWindow : Window
{
    private const int W = 512, H = 512;
    private const int SamplesPerPixel = 100;
    private readonly WriteableBitmap _wb;

    /// <summary>
    /// Initializes a new instance of the MainWindow class.
    /// Sets up the rendering surface and starts the path tracing process.
    /// </summary>
    public MainWindow()
    {
        InitializeComponent();
        _wb = new WriteableBitmap(W, H, 96, 96, PixelFormats.Bgra32, null);
        Canvas.Source = _wb;
        RenderScene();
    }

    /// <summary>
    /// Builds the camera coordinate system and viewport dimensions.
    /// </summary>
    /// <param name="E">Camera position (eye point).</param>
    /// <param name="lookAt">Point the camera is looking at.</param>
    /// <param name="fovDegrees">Vertical field of view in degrees.</param>
    /// <param name="aspect">Aspect ratio (width/height).</param>
    /// <returns>Camera basis vectors and viewport dimensions.</returns>
    /// <remarks>
    /// Returns a tuple containing:
    /// - f: forward direction (camera z-axis)
    /// - r: right direction (camera x-axis)
    /// - u: up direction (camera y-axis)
    /// - vw: viewport width
    /// - vh: viewport height
    /// </remarks>
    private static (Vector3 f, Vector3 r, Vector3 u, float vw, float vh) 
        BuildCameraBasis(Vector3 E, Vector3 lookAt, float fovDegrees, float aspect)
    {
        Vector3 f = Vector3.Normalize(lookAt - E);
        Vector3 upWorld = new Vector3(0, 1, 0);
        
        // Handle case where forward is parallel to world up
        if (MathF.Abs(Vector3.Dot(f, upWorld)) > 0.999f) 
            upWorld = new Vector3(0, 0, 1);
        
        Vector3 r = Vector3.Normalize(Vector3.Cross(upWorld, f));
        Vector3 u = Vector3.Cross(f, r);

        float fovRad = MathF.PI * fovDegrees / 180.0f;
        float vh = 2.0f * MathF.Tan(fovRad * 0.5f);
        float vw = vh * aspect;

        return (f, r, u, vw, vh);
    }

    /// <summary>
    /// Creates a camera ray for a specific pixel.
    /// </summary>
    /// <param name="x">Pixel x-coordinate.</param>
    /// <param name="y">Pixel y-coordinate.</param>
    /// <param name="W">Image width.</param>
    /// <param name="H">Image height.</param>
    /// <param name="E">Camera position.</param>
    /// <param name="cam">Camera basis and viewport information.</param>
    /// <returns>A ray from the camera through the specified pixel.</returns>
    private static Ray CreateEyeRay(int x, int y, int W, int H,
        Vector3 E, (Vector3 f, Vector3 r, Vector3 u, float vw, float vh) cam)
    {
        // Convert pixel coordinates to normalized device coordinates [-1, 1]
        float xNorm = (x / (float)(W - 1)) * 2.0f - 1.0f;
        float yNorm = ((H - 1 - y) / (float)(H - 1)) * 2.0f - 1.0f;

        float xOffset = xNorm * (cam.vw * 0.5f);
        float yOffset = yNorm * (cam.vh * 0.5f);

        Vector3 d = cam.f + xOffset * cam.r + yOffset * cam.u;
        d = Vector3.Normalize(d);

        return new Ray { E = E, d = d };
    }

    /// <summary>
    /// Renders the textured scene.
    /// </summary>
    private void RenderScene()
    {
        var scene = new Scene();
        
        string basePath = System.IO.Path.Combine(AppDomain.CurrentDomain.BaseDirectory, @"..\..\..\");
        
        Texture texture1 = null;
        Texture texture2 = null;
        Texture texture3 = null;
        Texture texture4 = null;
        Texture texture5 = null;
        Texture texture6 = null;
        Texture texture7 = null;
        Texture texture8 = null;
        Texture texture9 = null;
        Texture texture10 = null;
        Texture texture11 = null;
        Texture texture12 = null;
        
        texture1 = new Texture(System.IO.Path.Combine(basePath, "texture1.jpg"));
        texture2 = new Texture(System.IO.Path.Combine(basePath, "texture2.jpg"));
        texture3 = new Texture(System.IO.Path.Combine(basePath, "texture3.jpg"));
        texture4 = new Texture(System.IO.Path.Combine(basePath, "texture4.jpg"));
        texture5 = new Texture(System.IO.Path.Combine(basePath, "texture5.jpg"));
        texture6 = new Texture(System.IO.Path.Combine(basePath, "texture6.jpg"));
        texture7 = new Texture(System.IO.Path.Combine(basePath, "texture7.jpg"));
        texture8 = new Texture(System.IO.Path.Combine(basePath, "texture8.jpg"));
        texture9 = new Texture(System.IO.Path.Combine(basePath, "texture9.jpg"));
        texture10 = new Texture(System.IO.Path.Combine(basePath, "texture10.jpg"));
        texture11 = new Texture(System.IO.Path.Combine(basePath, "texture11.jpg"));
        texture12 = new Texture(System.IO.Path.Combine(basePath, "texture12.jpg"));
        

        // ========== SCENE SETUP ==========

        // Large emissive ceiling (area light) - positioned high above scene
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(0, 1010f, 0),  // Center much higher
            r = 1000f,                      // Large radius
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f), 
                Emission = new Vector3(2f, 2f, 2f),  // Bright white area light
                Roughness = 1.0f,
                Texture = null
            }
        });

        // Ground plane (large sphere below)
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(0, -1000.5f, 0),  // Center much lower
            r = 1000f,                         // Large radius
            Material = new Material 
            { 
                Color = new Vector3(0.3f, 0.3f, 0.3f),  // Medium gray to show shadows
                Emission = Vector3.Zero,
                Roughness = 1.0f,
                Texture = null
            }
        });
        
        // Sphere 0
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(12f, 0f, 2f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 8f, 
            Material = new Material 
            { 
                Color = new Vector3(.8f, .8f, .8f),  // White tint (multiplied with texture)
                Emission = Vector3.Zero,
                Roughness = 0.0f,  // Diffuse
                Texture = null
            }
        });

        // ========== TEXTURED SPHERES ==========

        // Sphere 1: Diffuse textured
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(-1.5f, 0f, 2f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 0.5f, 
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f),  // White tint (multiplied with texture)
                Emission = Vector3.Zero,
                Roughness = 1.0f,  // Diffuse
                Texture = texture1
            }
        });

        // Sphere 2: Glossy textured
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(0f, 0f, 2f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 0.5f, 
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f),
                Emission = Vector3.Zero,
                Roughness = 1.0f,  // Diffuse
                Texture = texture2
            }
        });

        // Sphere 3: Mirror-like textured
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(1.5f, 0f, 2f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 0.5f, 
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f),
                Emission = Vector3.Zero,
                Roughness = 1.0f,  // Diffuse
                Texture = texture3
            }
        });
        
        // Sphere 4: Pizza textured
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(-0.75f, 0f, 3.5f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 0.5f, 
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f),
                Emission = Vector3.Zero,
                Roughness = 1.0f,  // Diffuse
                Texture = texture4
            }
        });
        
        // Sphere 5
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(0.75f, 0f, 3.5f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 0.5f, 
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f),
                Emission = Vector3.Zero,
                Roughness = 1.0f,  // Diffuse
                Texture = texture5
            }
        });
        
        // Sphere 6
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(2.25f, 0f, 3.5f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 0.5f, 
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f),
                Emission = Vector3.Zero,
                Roughness = 1.0f,  // Diffuse
                Texture = texture6
            }
        });
        
        // Sphere 7
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(0f, 0f, 5f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 0.5f, 
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f),
                Emission = Vector3.Zero,
                Roughness = 1.0f,  // Diffuse
                Texture = texture7
            }
        });
        
        // Sphere 8
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(1.5f, 0f, 5f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 0.5f, 
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f),
                Emission = Vector3.Zero,
                Roughness = 1.0f,  // Diffuse
                Texture = texture8
            }
        });
        
        // Sphere 9
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(3f, 0f, 5f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 0.5f, 
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f),
                Emission = Vector3.Zero,
                Roughness = 1.0f,  // Diffuse
                Texture = texture9
            }
        });
        
        // Sphere 10
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(0.75f, 0f, 6.5f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 0.5f, 
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f),
                Emission = Vector3.Zero,
                Roughness = 1.0f,  // Diffuse
                Texture = texture10
            }
        });
        
        // Sphere 11
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(2.25f, 0f, 6.5f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 0.5f, 
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f),
                Emission = Vector3.Zero,
                Roughness = 1.0f,  // Diffuse
                Texture = texture11
            }
        });
        
        // Sphere 12
        scene.Spheres.Add(new Sphere 
        { 
            C = new Vector3(3.75f, 0f, 6.5f),  // Center at y=0, bottom at y=-0.5 (touching ground)
            r = 0.5f, 
            Material = new Material 
            { 
                Color = new Vector3(1f, 1f, 1f),
                Emission = Vector3.Zero,
                Roughness = 1.0f,  // Diffuse
                Texture = texture12
            }
        });
        

        // Camera setup - diagonal view, no ceiling visible
        Vector3 E = new Vector3(-2.5f, 2.8f, -0.5f);  // Diagonal position (left and front)
        Vector3 lookAt = new Vector3(0f, 0f, 2f);     // Look at center sphere
        float fov = 60f;
        
        float aspect = (float)W / H;
        var cam = BuildCameraBasis(E, lookAt, fov, aspect);

        var tracer = new PathTracer();
        int stride = _wb.BackBufferStride;
        byte[] pixels = new byte[stride * H];

        // Render with Monte Carlo sampling
        for (int y = 0; y < H; y++)
        {
            for (int x = 0; x < W; x++)
            {
                Vector3 colorAccum = Vector3.Zero;
                
                for (int s = 0; s < SamplesPerPixel; s++)
                {
                    Ray ray = CreateEyeRay(x, y, W, H, E, cam);
                    colorAccum += tracer.Trace(ray, scene);
                }

                Vector3 color = colorAccum / SamplesPerPixel;

                int idx = y * stride + x * 4;
                WriteBGRA(pixels, idx, color);
            }

        }

        _wb.Lock();
        _wb.WritePixels(new Int32Rect(0, 0, W, H), pixels, stride, 0);
        _wb.AddDirtyRect(new Int32Rect(0, 0, W, H));
        _wb.Unlock();

    }

    /// <summary>
    /// Converts a linear RGB color component to an 8-bit sRGB value.
    /// </summary>
    /// <param name="linear">Linear RGB component value (typically in range [0, 1]).</param>
    /// <returns>8-bit sRGB value in range [0, 255].</returns>
    /// <remarks>
    /// The sRGB transfer function applies a piecewise transformation:
    /// - For small values (≤ 0.0031308): linear scaling by 12.92
    /// - For larger values: power function with gamma ≈ 2.4
    /// This matches the standard sRGB specification (IEC 61966-2-1).
    /// </remarks>
    private static byte ToSrgb8(float linear)
    {
        if (linear <= 0f) return 0;
        if (linear >= 1f) return 255;

        float srgb = linear <= 0.0031308f
            ? 12.92f * linear
            : 1.055f * MathF.Pow(linear, 1f / 2.4f) - 0.055f;

        return (byte)Math.Clamp((int)MathF.Round(srgb * 255f), 0, 255);
    }
    
    /// <summary>
    /// Writes a linear RGB color to a byte buffer in BGRA32 format.
    /// </summary>
    /// <param name="buf">The byte array to write to (pixel buffer).</param>
    /// <param name="idx">The starting index in the buffer (must be 4-byte aligned).</param>
    /// <param name="linearRgb">The color in linear RGB space (components typically in [0, 1]).</param>
    /// <remarks>
    /// BGRA32 format stores pixels as: [Blue, Green, Red, Alpha] bytes.
    /// This is the native format for many graphics APIs and provides optimal performance.
    /// </remarks>
    private static void WriteBGRA(byte[] buf, int idx, Vector3 linearRgb)
    {
        buf[idx + 0] = ToSrgb8(linearRgb.Z); // Blue
        buf[idx + 1] = ToSrgb8(linearRgb.Y); // Green
        buf[idx + 2] = ToSrgb8(linearRgb.X); // Red
        buf[idx + 3] = 255;                   // Alpha
    }
}
