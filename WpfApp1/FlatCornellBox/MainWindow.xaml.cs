using System;
using System.Collections.Generic;
using System.Numerics;
using System.Windows;
using System.Windows.Media;
using System.Windows.Media.Imaging;

namespace FlatCornellBox;

/// <summary>
/// Represents a ray in 3D space for path tracing.
/// </summary>
public struct Ray
{
    public Vector3 E; // Origin
    public Vector3 d; // Direction
}

/// <summary>
/// Defines the material properties of a surface.
/// For flat shading, only the color property is used.
/// </summary>
public struct Material
{
    /// <summary>
    /// The albedo/diffuse color in linear RGB space.
    /// </summary>
    public Vector3 Color;
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
    public Material Material;
}

/// <summary>
/// Represents a sphere in the 3D scene.
/// Spheres are used both for actual spherical objects and as large primitives to create walls.
/// </summary>
public class Sphere
{
    public Vector3 C; // Center
    public float r; // Radius
    public Material Material;

    /// <summary>
    /// Computes the intersection of a ray with this sphere.
    /// Uses the quadratic formula to solve the ray-sphere equation.
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

        return new HitInfo
        {
            DidHit = true,
            distance = tHit,
            H = H,
            n = n,
            Material = this.Material
        };
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
/// Main window for the Flat Cornell Box renderer.
/// Implements a simple ray tracer that displays material colors without lighting calculations.
/// </summary>
public partial class MainWindow : Window
{
    private const int W = 512, H = 512;
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
    /// Creates the Cornell box scene and renders it using ray tracing.
    /// The Cornell box is constructed using large spheres as walls.
    /// </summary>
    /// <remarks>
    /// The scene consists of:
    /// - Red left wall
    /// - Blue right wall
    /// - Gray, floor, and back wall
    /// - White ceiling
    /// - Yellow sphere 
    /// - Cyan sphere
    /// </remarks>
    private void RenderScene()
    {
        // Create the scene
        var scene = new Scene();

        // Cornell box walls using large spheres
        // Left wall (red)
        scene.Spheres.Add(new Sphere
        {
            C = new Vector3(-1001, 0, 0),
            r = 1000,
            Material = new Material { Color = new Vector3(1f, 0f, 0f) }
        });

        // Right wall (blue)
        scene.Spheres.Add(new Sphere
        {
            C = new Vector3(1001, 0, 0),
            r = 1000,
            Material = new Material { Color = new Vector3(0f, 0f, 1f) }
        });

        // Ceiling (gray)
        scene.Spheres.Add(new Sphere
        {
            C = new Vector3(0, 1001, 0),
            r = 1000,
            Material = new Material { Color = new Vector3(1f, 1f, 1f) }
        });

        // Floor (gray)
        scene.Spheres.Add(new Sphere
        {
            C = new Vector3(0, -1001, 0),
            r = 1000,
            Material = new Material { Color = new Vector3(0.6f, 0.6f, 0.6f) }
        });

        // Back wall (gray)
        scene.Spheres.Add(new Sphere
        {
            C = new Vector3(0, 0, 1001),
            r = 1000,
            Material = new Material { Color = new Vector3(0.6f, 0.6f, 0.6f) }
        });

        // Yellow sphere
        scene.Spheres.Add(new Sphere
        {
            C = new Vector3(-0.6f, -0.7f, -0.6f),
            r = 0.3f,
            Material = new Material { Color = new Vector3(1f, 1f, 0f) }
        });

        // Cyan sphere
        scene.Spheres.Add(new Sphere
        {
            C = new Vector3(0.3f, -0.4f, 0.3f),
            r = 0.6f,
            Material = new Material { Color = new Vector3(0.8f, 0.9f,0.9f) }
        });
        
        // ==================================================
        // ================== Camera setup ==================
        // ==================================================
        
        // Camera setup - Perspective 1
        /*Vector3 eye = new Vector3(0f, 0f, -4.5f);
        float fov = 36.0f;
        Vector3 lookAt = new Vector3(0, 0, 6);*/

        // Camera setup - Perspective 2
        Vector3 eye = new Vector3(-0.9f, -0.5f, 0.9f);
        float fov = 110.0f;
        Vector3 lookAt = new Vector3(0, 0, 0);
        
        Vector3 backgroundColor = new Vector3(0, 0, 0);

        // Compute camera basis vectors
        float aspectRatio = (float)W /H;
        float fovRad = MathF.PI /180.0f * fov;
        float viewportHeight = 2.0f * MathF.Tan(fovRad / 2.0f);
        float viewportWidth = viewportHeight * aspectRatio;

        Vector3 forward = Vector3.Normalize(lookAt - eye);
        Vector3 right = Vector3.Normalize(Vector3.Cross(new Vector3(0, 1, 0), forward));
        Vector3 up = Vector3.Cross(forward, right);

        // Prepare pixel buffer
        int stride = _wb.BackBufferStride;
        byte[] pixels = new byte[stride *H];

        // Render each pixel
        for (int y = 0;y < H; y++)
        {
            for (int x = 0; x < W; x++)
            {
                // Convert pixel coordinates to normalized device coordinates [-1, 1]
                float u = (x / (float)(W - 1)) * 2.0f - 1.0f;
                float v = ((H - 1 - y)/ (float)(H - 1)) * 2.0f - 1.0f;

                // Create eye ray for this pixel
                Ray ray = new Ray
                {
                    E = eye,
                    d = Vector3.Normalize(
                        forward + u * viewportWidth /2 *right + v * viewportHeight /2 * up)
                };

                // Find closest intersection
                HitInfo hit = scene.Intersect(ray);

                // Compute color (flat shading: just return material color)
                Vector3 linearColor = hit.DidHit ?hit.Material.Color :backgroundColor;

                // Write pixel to buffer
                int idx = y * stride + x *4;
                WriteBGRA(pixels, idx, linearColor);
            }
        }

        // Update the bitmap
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